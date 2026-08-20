"""
======================COPYRIGHT/LICENSE START==========================

Validation.py: Utility functions to work with validation objects. 

Copyright (C) 2004 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue. The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Accepted by Proteins (2004).

===========================REFERENCE END===============================

"""

VALID_VALUE_TYPE_ATTRS = {type(0.0): 'floatValue',
                          type(1): 'intValue',
                          type('a'): 'textValue',
                          type(True): 'booleanValue'} 


from ccp.util.Software import getSoftware

# # # # # # #  C I N G  E X A M P L E  # # # # # # #

def storeRogScores(ccpnEnsemble, scores, context='CING', cingVersion='1.0'):
   # Assumes scores in same order as residues
   
   keyword    = 'ROG'
   definition = 'Overall per-residue validation score for display. ROG=Red/Orange/Green'
   synonym    = 'Residue ROG Score'
   
   project  = ccpnEnsemble.root
   software = getSoftware(project, context, cingVersion)
   
   validStore = getEnsembleValidationStore(ccpnEnsemble, context,
                                           keywords=[keyword, ],
                                           definitions=[definition, ],
                                           synonyms=[synonym, ],
                                           software=software)
                                           
   residues = []
   for chain in ccpnEnsemble.sortedCoordChains():
     residues.extend(chain.sortedResidues())
   
   storeResidueValidations(validStore, context, keyword, residues, scores)

# # # # # # # # # # # # # # # # # # # # # # # # # # #



def getEnsembleValidationStore(ensemble, context, keywords=None,
                               definitions=None, synonyms=None,
                               software=None):
  """Descrn: Get a CCPN object to store validation results for an ensemble
             in a given program context. Requires a list of keywords which will
             be used in this context. Allows optional lists of definitions and
             user-friendly synonyms for these keywords.
             Optional argument for passing software specification (otherwise
             defaults to the name of the context)
     Inputs: MolStructure.StructureEnsemble, Word, List of Words,
             List of Lines, List of Words, Method.Software  
     Output: Validation.ValidationStore
  """
  
  if not keywords:
    keywords = []

  memopsRoot = ensemble.root
  eid = '%s_%s' % (context, ensemble.guid)
  validStore = ensemble.findFirstValidationStore(name=eid)
  
  if validStore is None:
    if not software:
      software = getSoftware(memopsRoot, context, '?')
    
    validStore = memopsRoot.newValidationStore(name=eid, software=software,
                                               structureEnsemble=ensemble)
  
  validStore.nmrProject = memopsRoot.currentNmrProject
  
  keywordStore = memopsRoot.findFirstKeywordDefinitionStore(context=context)
  
  if not keywordStore:
    keywordStore = memopsRoot.newKeywordDefinitionStore(context=context)
  
  for i, keyword in enumerate(keywords):
    keywordDefinition = keywordStore.findFirstKeywordDefinition(keyword=keyword)
    
    if not keywordDefinition:
      keywordDefinition = keywordStore.newKeywordDefinition(keyword=keyword)
    
    if definitions and (i < len(definitions)):
      keywordDefinition.explanation = definitions[i]
    
    if synonyms and (i < len(synonyms)):
      keywordDefinition.name = synonyms[i]
          
  return validStore


def getValidationObjects(validStore, className, context, keyword):
  """Descrn: Find a given class of validation objects in a validation store
             in a given (program) context, with a given keyword.
     Inputs: Validation.ValidationStore, Word, Word, Word
     Output: List of Validation.Validations
  """
  return validStore.findAllValidationResults(context=context,
                                             keyword=keyword, 
                                             className=className)


def replaceValidationObjects(validStore, className, keyword, context, dataList):
  """Descrn: Store validation data as CCPN validation objects,
             overwriting all previous records of such information
             in the validation store object. Finds a given class of
             validation objects a given keyword in a given context.
             The input data list is a 2-tuple containing a list of the
             CCPN objects validated and the value associated with that
             validation. Note that the validated CCPN object types must
             match the type required by the validation object className.
             *NOTE* this is a slow function and should often be replaced
             with a class-specific equivalent.
     Inputs: Validation.ValidationStore, Word, Word, Word,
             List of 2-Tuples of (List of CCPN objects, )
     Output: None
  """

  for validObj in getValidationObjects(validStore, className, keyword, context):
    validObj.delete()

  newObject = getattr(validStore, 'new%s' % className)

  validatedObjAttr = None

  for validatedObjects, value in dataList:
    validObj = newObject(context=context, keyword=keyword)
    
    if not validatedObjAttr:
      for role in validObj.metaclass.getAllRoles():
        if role.locard == 0:
          validatedObjAttr = role.name
    
    if validatedObjAttr:
      setattr(validObj, validatedObjAttr, validatedObjects)
    
    valueAttr = VALID_VALUE_TYPE_ATTRS.get(type(value), 'textValue')
    setattr(validObj, valueAttr, value)
    
  
def getResidueValidation(validStore, residue, context, keyword):
  """Descrn: Get any existing residue validation results from a CCPN
             validation store which have the given keywords in the
             given (program) context.
             *NOTE* This function may be quicker than using the generic
             getValidationObjects() because the link is queried from the
             validated object, not the validation store, which often
             has fewer total validation objects.
     Inputs: Validation.ValidationStore,
             MolStructure.Residue, Word, Word
     Output: Validation.ResidueValidation  
  """
  
  # Define data model call to find exting result
  findValidation = residue.findFirstResidueValidation
  
  validObj = findValidation(validationStore=validStore,
                            context=context, keyword=keyword)

  return validObj
 

def storeResidueValidations(validStore, context, keyword, residues, scores):
  """Descrn: Store the per-residue scores for a an ensemble within
             CCPN validation objects.
             *NOTE* This function may be quicker than using the generic
             replaceValidationObjects() because it is class specifc
     Inputs: Validation.ValidationStore,
             List of MolStructure.Residues, List of Floats
     Output: List of Validation.ResidueValidations  
  """
  
  validObjs = []
    
  # Define data model call for new result
  newValidation = validStore.newResidueValidation
 
  for i, residue in enumerate(residues):
  
    score = scores[i]
    
    # Find any existing residue validation objects
    validObj = getResidueValidation(validStore, residue, context, keyword)

    # Validated object(s) must be in a list
    residueObjs = [residue, ]
 
    # Make a new validation object if none was found
    if not validObj:
      validObj = newValidation(context=context, keyword=keyword,
                               residues=residueObjs)

    # Set value of the score
    validObj.floatValue = score
    
    validObjs.append(validObj)
  
  return validObjs

def getModelValidation(validStore, model, context, keyword):
  """Descrn: Get any existing structure validation results from a CCPN
             validation store which have the given keywords in the
             given (program) context. Finds only structure validations
             that refer to a single model.
             *NOTE* This function may be quicker than using the generic
             getValidationObjects() because the link is queried from the
             validated object, not the validation store, which often
             has fewer total validation objects.
     Inputs: Validation.ValidationStore,
             MolStructure.Model, Word, Word
     Output: Validation.StructureValidation  
  """
  
  # Define data model call to find exting result
  findValidation = model.findFirstStructureValidation
  
  validObj = findValidation(validationStore=validStore, models=set([model,]),
                            context=context, keyword=keyword)

  return validObj

def storeModelValidations(validStore, context, keyword, models, scores):
  """Descrn: Store the per-model scores for a an ensemble within
             CCPN validation objects. Creates one validation object
             for each model separately (not a combination).
             *NOTE* This function may be quicker than using the generic
             replaceValidationObjects() because it is class specifc
     Inputs: Validation.ValidationStore, Word, Word,
             List of MolStructure.Models, List of Floats
     Output: List of Validation.StructureValidations  
  """
  
  validObjs = []
    
  # Define data model call for new result
  newValidation = validStore.newStructureValidation
 
  for i, model in enumerate(models):
  
    score = scores[i]
    
    # Find any existing validation objects for model
    validObj = getModelValidation(validStore, model, context, keyword)

    # Validated object must be in a list
    modelObjs = [model, ]
 
    # Make a new validation object if none was found
    if not validObj:
      validObj = newValidation(context=context, keyword=keyword,
                               models=modelObjs)

    # Set value of the score
    validObj.floatValue = score
    
    validObjs.append(validObj)
  
  return validObjs


def storeEnsembleValidation(validStore, context, keyword, ensemble, score):
  """Descrn: Store scores for an ensemble as CCPN validation objects.
     Inputs: Validation.ValidationStore, Word, Word,
             MolStructure.StructureEnsemble, Float
     Output: Validation.EnsembleValidation    
  """
  
  findValidation = ensemble.findFirstEnsembleValidation

  validObj = findValidation(validationStore=validStore,
                            context=context, keyword=keyword)
  
  if not validObj:
    validObj = validStore.newEnsembleValidation(context=context, keyword=keyword)
  
  validObj.floatValue = score
  
  return validObj

