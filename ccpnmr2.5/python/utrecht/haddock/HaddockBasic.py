#!/usr/bin/env python

"""
=========================================================================
Package:    - Code for the Graphical User Interface fontend to the 
              Haddock model package in the CCPN data model. 
            - Code for the export of a Haddock compatible project. A 
              Haddock compatible project can either be a parameter file
              ready for submission to the Haddock webserver or a
              directory structure with necessary files for use with a 
              localy installed version of Haddock.

Dependencies: The CCPN Haddock package requires CCPN data model version
              2.0 or higher. The export of a webserver compatible 
              parameter file requires Haddock webserver version 2.1 or 
              higher and a valid user account. The export of a 'classic' 
              Haddock project requires Haddock version 2.0 or higher.

Copyright and License information:
              The Haddock data model as implemented in the CCPN data
              model as well as the use of CCPN GUI code elements is 
              licenced to the CCPN Projects (Copyright (C) 2008) and
              distributed under the terms of the GNU Lesser General
              Public License.
            
              The Haddock project export code as well as the use of 
              Haddock software is covert in the Haddock License
              agreement (Copyright (C) 2008 Haddock Project, Bijvoet
              Center for Biomolecular Research, Utrecht University,
              The Netherlands).

GNU LGPL:        This library is free software; you can redistribute it 
              and/or modify it under the terms of the GNU Lesser General 
              Public License as published by the Free Software 
              Foundation; either version 2.1 of the License, or (at 
              your option) any later version.
 
              A copy of this license can be found in LGPL.license
 
              This library is distributed in the hope that it will be 
              useful, but WITHOUT ANY WARRANTY; without even the implied 
              warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
              PURPOSE. See the GNU Lesser General Public License for 
              more details.
 
              You should have received a copy of the GNU Lesser General 
              Public License along with this library; if not, write to 
              the Free Software Foundation, Inc., 59 Temple Place, Suite 
              330, Boston, MA 02111-1307 USA.

Information:  For further information regarding CCPN, please contact:
              - CCPN website (http://www.ccpn.ac.uk/)
              - email: ccpn@bioc.cam.ac.uk
              
              For further information regarding Haddock, please contact
              Alexandre M.J.J. Bonvin:
              - http://haddock.chem.uu.nl
              - email: a.m.j.j.bonvin@uu.nl    

Citing:          If you are using this software for academic purposes, we 
                suggest quoting the following references:

              For CCPN:    
              Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
              Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. 
              Ionides and Ernest D. Laue (2005). A framework for 
              scientific data modeling and automated software development. 
              Bioinformatics 21, 1678-1684.
            
              For Haddock:
              Cyril Dominguez, Rolf Boelens and Alexandre M.J.J. Bonvin 
              (2003). HADDOCK: a protein-protein docking approach based 
              on biochemical and/or biophysical information. 
              J. Am. Chem. Soc. 125, 1731-1737.
            
              S.J. de Vries, A.D.J. van Dijk, M. Krzeminski, M. van Dijk, 
              A. Thureau, V. Hsu, T. Wassenaar and A.M.J.J. Bonvin (2007) 
              HADDOCK versus HADDOCK: New features and performance of 
              HADDOCK2.0 on the CAPRI targets. 
              Proteins: Struc. Funct. & Bioinformatic 69, 726-733.    
=========================================================================
"""

#import sys
#import math
#from os.path     import join, exists, split, isdir, isfile
#from os          import rename

import os

#from ccp.general.Io                      import getStdChemComps, getChemComp
#from ccpnmr.analysis.core.MoleculeBasic  import findMatchingChain, getLinkedResidue
#from ccpnmr.analysis.core.StructureBasic import getBestNamingSystem, findMatchingMolSystemAtom
from ccpnmr.format.converters.CnsFormat  import CnsFormat
#from ccp.util.Molecule                   import makeChain, addMolResidues, makeMolecule, nextChainCode
#from memops.gui.MessageReporter          import showOkCancel,showWarning,showYesNo
#from memops.universal.Util               import returnInt, returnFloat

# New import, replaces previous code duplication
from ccp.lib.StructureIo import getStructureFromFile

from HaddockLocal                        import rdcProtocolStore, daniProtocolStore

def addRdcParam(run,termId):

    """Add an RDC energy protocol to the RDC energyTerm. Both termId's need to match so that the RDC file can
       be linked to its RDC energy protocol
    """
    energyTermStore = run.newHaddockEnergyTerm(code='rdcProtocolStore',termId=termId)
    
    terms = rdcProtocolStore['terms'].keys()
    for term in terms:
        energyTerm = energyTermStore.newEnergyTermParameter(code=term,value=rdcProtocolStore['terms'][term])

def addDaniParam(run,termId):

    """Add an DANI energy protocol to the DANI energyTerm. Both termId's need to match so that the DANI file can
       be linked to its DANI energy protocol
    """
    energyTermStore = run.newHaddockEnergyTerm(code='daniProtocolStore',termId=termId)

    terms = daniProtocolStore['terms'].keys()
    for term in terms:
        energyTerm = energyTermStore.newEnergyTermParameter(code=term,value=daniProtocolStore['terms'][term])

def getPdbString(model,chains,chainRename=None,blankchain=False,fileEnd='END'):

    """Description:    Convert a CCPN molEnsemble into a PDB formatted string Atom names will use IUPAC nomenclature
       Input:         MolStructure.StructureEnsemble, list of chains (pdb One Letter Code) to check if we need to 
                    write only a subset of chains of the MolStructure. The chainID can be set to the haddock partner
                    code by using the chainRename argument. If the chain argument is set to True both chain and segid
                    are printed (needed for parameter file).
       Output:         PDB structure String
    """
    project = model.root
    cnsFormat = CnsFormat(project,guiParent=None)

    if chainRename: forceExportSegId = chainRename
    else: forceExportSegId = " "

    if blankchain: forceExportChainId = " "
    else: forceExportChainId = forceExportSegId
    
    if not chains: chainlist = project.findFirstMolSystem().sortedChains()
    
    cnsFormat.writeCoordinates('void',
                                structures = [model],
                                forceExportSegId = forceExportSegId,
                                forceExportChainId = forceExportChainId,
                                resetMapping=True,
                                exportChains = chains,
                                minimalPrompts = True,
                                noWrite = True)

    cnsFormat.coordinateFile.write(endStatement=fileEnd,writeString=True)

    return cnsFormat.coordinateFile.coordFileString

def getAirSegments(haddockPartner):
    
    """Description: Get lists of the active and passive residues for ambiguous
                       interaction restraints from a haddock interacting partner.
                       Numbers output are in the HADDOCK numbering system for the
                       partner (i.e. may be different from CCPN chains)
       Input:         Haddock.HaddockPartner
       Output:         Dict of Lists (Dict['active'/'passive'] = List of residues)
    """
    segmentDict = {'active':[],'passive':[]}
    
    for chain in haddockPartner.chains:
        residues = [(r.residue.seqCode, r) for r in chain.residues]
        residues.sort()
        
        for seqId, residue in residues:
            interaction = residue.interaction
            if interaction == 'active': segmentDict['active'].append(seqId)
            elif interaction == 'passive': segmentDict['passive'].append(seqId)
    
    return segmentDict

def getFlexibleResidues(haddockPartner,export='zone'):
    
    """Description: Get lists of the flexible and semi-flexible residues for a
                       haddock interacting partner. Numbers output are in the HADDOCK
                       numbering system for the partner (i.e. may be different from CCPN chains)
       Arguments:    export = 'zone', will return semi and fully flexible residues as a list
                    of tuples containing start and end residues of residue zones.
                    export = 'list', will return semi and fully flexible residues as a plain
                    list of residues. 
       Input:         Haddock.HaddockPartner
       Output:         Dict of Lists (Dict['semi'/'full'] = List of residues)
    """
    segmentDict = {'semi':[],'full':[]}
    
    for chain in haddockPartner.chains:
        residues = [(r.residue.seqCode, r) for r in chain.residues]
        residues.sort()
        
        for seqId, residue in residues:
            flexibility = residue.flexibility
            if flexibility == 'semi': segmentDict['semi'].append(seqId)
            elif flexibility == 'full': segmentDict['full'].append(seqId)
    
    if export == 'zone':
        for flex in segmentDict:
            if len(segmentDict[flex]):
                flexzones = []
                flexresidues = segmentDict[flex]
                flexresidues.sort()
                start = flexresidues[0]
                for r in range(len(flexresidues)):
                    if not r == len(flexresidues)-1:
                        if not flexresidues[r+1] == flexresidues[r]+1:
                            flexzones.append((start,flexresidues[r]))
                            start = flexresidues[r+1]
                    else:
                        flexzones.append((start,flexresidues[r]))
                segmentDict[flex] = flexzones
        return segmentDict         
    else: return segmentDict

def makeBackup(infile):
    
    """Make backup of files or directories by checking if file (or backup as _*) is there and rename"""
    
    if os.path.isfile(infile) or os.path.isdir(infile):
        i = 1
        while i < 100:
            if os.path.isfile(infile+'_'+str(i)) or os.path.isdir(infile+'_'+str(i)): i += 1
            else:
                os.rename(infile, infile+'_'+str(i))
                print("File %s exists\nrename to %s" % (infile, infile+'_'+str(i)))
                break

def copyRun(run, nmrConstraintStore=None):
    """Descrn: Make a HADDOCK run based upon an existing one. Copies existing
               attributes, links and children (makes equivalent)
       Inputs: Hadock.Run, NmrConstraint.NmrConstraintStore
       Output: Hadock.Run
    """

    project = run.haddockProject

    runB = project.newRun(nmrConstraintStore=nmrConstraintStore)

    for metaAttr in run.metaclass.attributes:
        if metaAttr.changeability == 'changeable':
            attrName = metaAttr.name
            value = getattr(run,attrName)
            setattr(runB,attrName,value)

    for symmetry in run.symmetryRestraints:
        runB.addSymmetryRestraint(symmetry)

    for scoringWeight in run.scoringWeights:
        runB.newScoringWeight(stage=scoringWeight.stage,
                                term=scoringWeight.term,
                                value=scoringWeight.value)

    for term in run.haddockEnergyTerms:
        termB = runB.newHaddockEnergyTerm(code=term.code,
                                        termId=term.termId,
                                        name=term.name,
                                        fileName=term.fileName,
                                        details=term.details)
        
        if run.nmrConstraintStore is nmrConstraintStore:
            termB.constraintList = term.constraintList
    
        for param in term.energyTermParameters:
            termB.newEnergyParameter(code=param.code, value=param.value)

    return runB

def setRunConstraintSet(run, nmrConstraintStore):
    """Descrn: Set the constraint set for a haddock run. Sets a new one
               if none exists and otherwise if changing constraint sets
               the HaddockEnergy termns are checked and moved over for the
               different constraint lists.
       Inputs: Hadock.Run, NmrConstraint.NmrConstraintStore
       Output: None
    """

    if nmrConstraintStore is not run.nmrConstraintStore:
        for term in run.haddockEnergyTerms:
            term.constraintList = None
        
        run.__dict__['nmrConstraintStore'] = nmrConstraintStore
        
    return run

def setPartnerEnsemble(haddockPartner, ensemble):
    """Descrn: Set the structure ensemble for a haddock partner.
               Sets molSystem and chains as required.
       Inputs: Hadock.HaddockPartner, MolStructure.StructureEnsemble
       Output: None
    """

    haddockPartner.molSystem = ensemble.molSystem

    chains = [c.chain for c in ensemble.coordChains]

    setPartnerChains(haddockPartner, chains)

def setPartnerChains(haddockPartner, chains):
    """Descrn: Set the haddock chains (mappings to ccp molSystem chains) for a haddock partner
               using a list of MolSystem.Chains - adds and deleted haddock chains as needed.
               Also sets the haddock resiues appropriately.
       Inputs: Hadock.HaddockPartner, Lost of MolSystem.Chains
       Output: None
    """

    for hChain in haddockPartner.chains:
        if hChain.chain not in chains: hChain.delete()

    molType = None
    for chain in chains:
        if molType is None: molType = chain.molecule.molType
        elif molType != chain.molecule.molType:
            print 'CCPN-HADDOCK setPartnerChains failed: Chains not of same type'
            return

    for chain in chains:
        hChain = haddockPartner.findFirstChain(chain=chain)
        if not hChain: haddockPartner.newChain(chain=chain)

        molType = chain.molecule.molType

    if molType == 'DNA':
        haddockPartner.isDna = True
        haddockPartner.forceFieldCode = 'DNA'

    elif molType == 'RNA' :
        haddockPartner.isDna = True
        haddockPartner.forceFieldCode = 'RNA'

    else:
        haddockPartner.isDna = False
        haddockPartner.forceFieldCode = 'TOPALLHDG'

    # Curate all haddock residue numbers and CCPN residue links
    hSeqId = 1
    for hChain in haddockPartner.sortedChains():
        chain     = hChain.chain
        residues  = chain.sortedResidues()
        hResidues = hChain.sortedResidues()

        resDict = {}
        for hResidue in hResidues:
            residue = hResidue.residue
            resDict[residue] = hResidue

        for residue in residues:
            hResidue = resDict.get(residue)

            if hResidue is None: hResidue = hChain.newResidue(haddockSeqId=hSeqId,residue=residue)
            else: hResidue.haddockSeqId = hSeqId

            hSeqId += 1


class evalWcPairing(object):

    """Description: Class to evaluate Watson-Crick hydrogen bonding between
                    bases in regular RNA or DNA structures. WC-pairs are
                    evaluated between possible hydrogen bond donors and 
                    acceptors in bases. The nucleotide type does not matter.
                    An heavy atom upper distance limit of 3.0 A by default is
                    used as cutoff.
       Input      : ccpn.molecule.MolSystem instance
       Output      : List of tuples with ccpn.molecule.MolStructure.Residue
                    instances involved in WC pairing
       Arguments  : Hydrogen bond upper distance limit. 3.0 A by default
    """
    def __init__(self,partner=None,hbond=3.0):

        self.partner = partner.structureEnsemble
        self.allowed = [('O6','N4'),
                        ('N4','O6'),
                        ('N1','N3'),
                        ('N3','N1'),
                        ('N2','O2'),
                        ('O2','N2'),
                        ('O4','N6'),
                        ('N6','O4')]
        self.hbond = hbond**2    
        self.pairs = []                

        if partner.isDna:
            self.__makeAtomSelection__()
            self.__createDistanceMatrix__()
            self.__resolveWcPairs__()

    def __makeAtomSelection__(self):

        """For every residue in the system select the possible base hydrogen 
           bond donors and acceptors
        """
        self.atoms = []
        allowedAtoms = [n[0] for n in self.allowed]

        hChains = self.partner.sortedCoordChains()
        for hChain in hChains:
            for hResidue in hChain.sortedResidues():
                self.atoms += ([a for a in hResidue.sortedAtoms() if a.name in allowedAtoms])

    def __createDistanceMatrix__(self):

        """Calculate the distance between all possible base hydrogen bond donors and acceptors. Checks for:
           no distance calculation within same objects and between same base types"""

        self.matrix = []

        for atm1 in self.atoms:
            for atm2 in self.atoms:
                if not atm1.residue.residue.ccpCode == atm2.residue.residue.ccpCode and (atm1.name,atm2.name) in self.allowed:
                    dist = self.__evalAtomAtomDistance__(atm1.findFirstCoord(),atm2.findFirstCoord())
                    self.matrix.append((atm1,atm2,dist))

    def __evalAtomAtomDistance__(self,atm1,atm2):

        """Calculate the distance between atom1 and atom2. For speed we do not take the square root to optain
           the final distance.
        """
        cor1 = (atm1.x,atm1.y,atm1.z); cor2 = (atm2.x,atm2.y,atm2.z)

        d= float(0)
        if cor1 == None or cor2 == None: return None
        else:
            for i in range(len(cor1)): d=d+(cor1[i]-cor2[i])**2
            return d

    def __resolveWcPairs__(self):

        """Compose the final list of base-pairs based on the WC pairing profile."""

        for dist in self.matrix:
            if dist[2] <= self.hbond:
                pair1 = (dist[0].residue,dist[1].residue)
                pair2 = (dist[1].residue,dist[0].residue)
                if not pair1 in self.pairs and not pair2 in self.pairs:
                    self.pairs.append(pair1)
