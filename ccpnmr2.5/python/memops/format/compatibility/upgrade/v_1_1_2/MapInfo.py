""" NB this file is deliberately empty of real information
"""
# Packages, classElements and AbstractDataTypes skipped in new model
# (prefix, typeName, elemName, newGuid, elemType)
skipElements = [
]

# classElements skipped in new model, but available for simple data transfer
# (prefix, typeName, elemName, newGuid, elemMap, valueTypeGuid)
delayElements = [
]

# MetaConstraints added in new model
# (qualifiedName, guid)
newConstraints = [
]

# Mandatory classElements added in new model
# New ClassElements with locard !=0, no default, not derived or Implementation
# (prefix, typeName, elemName, newGuid)
newMandatories = [
 ('NMRX', 'SystematicName', 'name', 'www.ccpn.ac.uk_Fogh_2008-03-30-16:37:37_00004'), 
 ('NMRX', 'SystematicName', 'namingSystem', 'www.ccpn.ac.uk_Fogh_2008-03-30-16:37:37_00003'), 
 ('NMRX', 'SystematicName', 'refExperiment', 'www.ccpn.ac.uk_Fogh_2008-03-30-16:37:37_00001'), 
]

# Packages, classElements and AbstractDataTypes added in new model
# Optional, i.e. excluding mandatory classElements given above
# (prefix, typeName, elemName, newGuid)
newElements = [
 ('NMRX', 'NamingSystem', None, 'www.ccpn.ac.uk_Fogh_2008-03-30-16:37:31_00002'), 
 ('NMRX', 'RefExperiment', 'systematicNames', 'www.ccpn.ac.uk_Fogh_2008-03-30-16:37:37_00002'), 
 ('NMRX', 'SystematicName', None, 'www.ccpn.ac.uk_Fogh_2008-03-30-16:37:31_00001'), 
]

# Class elements that exist in both models but that require handcode for
# transfer. E.g. elements that go from derived to non-derived.
# Note that old derivation functions can not be relied on to work during
# data transfer
# (prefix, typeName, elemName, newGuid, elemType)
neutraliseElements = [
]

# Differences between equivalent classElements and AbstractDataTypes :

# name changes
# (prefix, typeName, elemName, newName, newGuid
renames = [
 ('DANG', 'DangleResidue', 'philower', 'phiLower', 'www.ccpn.ac.uk_Fogh_2007-10-01-14:14:02_00007'), 
 ('CHEM', 'ChemComp', 'glycoCtCode', 'baseGlycoCtCode', 'www.ccpn.ac.uk_Fogh_2007-04-02-13:31:25_00001'), 
 ('CHEM', 'ChemCompVar', 'varGlycoCtCode', 'glycoCtCode', 'www.ccpn.ac.uk_Fogh_2007-04-02-13:31:25_00003'), 
 ('NMRX', 'RefExpDimRef', 'groupingNo', 'groupingId', 'www.ccpn.ac.uk_Fogh_2007-04-03-15:13:51_00001'), 
]

# ValueType changes
# change types are : 'ignore': do nothing, 'delay': available for calculation
#                    'skip' not read.
# (prefix, typeName, elemName, action, newGuid, elemMap, valueTypeGuid)
typeChanges = [
]

# Different elements with matching qualifiedNames
# (element.qName, differentTags
nameMatches = [
 ('ccp.molecule.ChemComp.ChemCompVar.glycoCtCode', set(['isDerived', 'guid', 'documentation'])), 
]

# Differences for matching elements, 
# excluding those where only names and/or valueTypes differ
# (oldElem.qName, newElem.name, oldGuid, newGuid, differentTags
allDiffs = [
]
