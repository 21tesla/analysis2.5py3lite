
# Packages, classElements and AbstractDataTypes skipped in new model
# (prefix, typeName, elemName, newGuid, elemType)
skipElements = [
]

# classElements skipped in new model, but available for simple data transfer
# (prefix, typeName, elemName, newGuid, elemMap, valueTypeGuid)
delayElements = [
 ('ANAY', 'WindowPanelGroup', 'showSameSpectra', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:04:35_00004', {'tag': 'ANAY.WindowPanelGroup.showSameSpectra', 'type': 'attr', 'name': 'showSameSpectra'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00028'), 
 ('CALC', 'NmrCalcStore', 'methodStoreName', 'www.ccpn.ac.uk_Fogh_2010-05-10-13:46:58_00003', {'eType': 'cplx', 'tag': 'CALC.NmrCalcStore.methodStoreName', 'type': 'attr', 'name': 'methodStoreName'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00033'), 
 ('CALC', 'NmrCalcStore', 'softwareName', 'www.ccpn.ac.uk_Fogh_2010-05-10-13:46:58_00004', {'eType': 'cplx', 'tag': 'CALC.NmrCalcStore.softwareName', 'type': 'attr', 'name': 'softwareName'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00033'), 
 ('CALC', 'NmrCalcStore', 'softwareVersion', 'www.ccpn.ac.uk_Fogh_2010-05-10-13:46:58_00005', {'eType': 'cplx', 'tag': 'CALC.NmrCalcStore.softwareVersion', 'type': 'attr', 'name': 'softwareVersion'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00033'), 
]

# MetaConstraints added in new model
# (qualifiedName, guid)
newConstraints = [
 ('cambridge.WmsProtocol.InterfaceParameter.hicard.hicard_consistent_with_ProtocolParameter_hicard', 'www.ccpn.ac.uk_Fogh_2013-10-11-09:59:51_00001'), 
 ('cambridge.WmsProtocol.InterfaceParameter.locard.locard_consistent_with_ProtocolParameter_locard', 'www.ccpn.ac.uk_Fogh_2013-10-11-09:59:51_00002'), 
 ('ccp.nmr.NmrCalc.Run.masterRun.derived_runs_cannot_be_nested', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00003'), 
]

# Mandatory classElements added in new model
# New ClassElements with locard !=0, no default, not derived or Implementation
# (prefix, typeName, elemName, newGuid)
newMandatories = [
 ('NMRS', 'ExperimentWeight', 'expCode', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00001'), 
 ('NMRS', 'ExperimentWeight', 'trialSet', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00005'), 
]

# Packages, classElements and AbstractDataTypes added in new model
# Optional, i.e. excluding mandatory classElements given above
# (prefix, typeName, elemName, newGuid)
newElements = [
 ('ANAL', 'SpectrumWindow', 'isZeroLineShown', 'www.ccpn.ac.uk_Fogh_2013-05-03-11:50:01_00001'), 
 ('ANAL', 'SpectrumWindowView', 'isContourLineVisible', 'www.ccpn.ac.uk_Fogh_2013-05-07-17:07:06_00001'), 
 ('ANAW', 'AbstractModule', 'details', 'www.ccpn.ac.uk_Fogh_2012-09-10-14:34:37_00001'), 
 ('ANAW', 'AbstractModule', 'keywords', 'www.ccpn.ac.uk_Fogh_2012-09-10-14:34:37_00002'), 
 ('ANAY', 'AbstractPanel', 'layoutArea', 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:29_00001'), 
 ('ANAY', 'Layout', 'rank', 'www.ccpn.ac.uk_Fogh_2012-09-10-14:34:35_00001'), 
 ('ANAY', 'LayoutArea', None, 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:27_00002'), 
 ('ANAY', 'SpectrumSharing', None, 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:27_00001'), 
 ('ANAY', 'WindowPanel1d', 'showPeakPickLevel', 'www.ccpn.ac.uk_Fogh_2012-09-17-10:38:43_00001'), 
 ('ANAY', 'WindowPanelGroup', 'spectrumSharing', 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:29_00002'), 
 ('CALC', 'Run', 'derivedRuns', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00001'), 
 ('CALC', 'Run', 'masterRun', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00002'), 
 ('CALC', 'Run', 'methodStoreName', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00004'), 
 ('CALC', 'Run', 'softwareName', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00005'), 
 ('CALC', 'Run', 'softwareVersion', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00006'), 
 ('CALC', 'Run', 'wmsProtocolName', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00007'), 
 ('NMR', 'Experiment', 'userExpCode', 'www.ccpn.ac.uk_Fogh_2012-07-23-11:52:27_00001'), 
 ('NMRC', 'FixedResonance', 'covalentlyBound', 'www.ccpn.ac.uk_Fogh_2012-06-25-14:41:56_00001'), 
 ('NMRS', 'ExperimentWeight', None, 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:44_00001'), 
 ('NMRS', 'ExperimentWeight', 'changeIsPositive', 'www.ccpn.ac.uk_Fogh_2012-07-09-10:23:18_00001'), 
 ('NMRS', 'ExperimentWeight', 'intensityScale', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00003'), 
 ('NMRS', 'ExperimentWeight', 'meritThreshold', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00002'), 
 ('NMRS', 'ExperimentWeight', 'weight', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00004'), 
 ('NMRS', 'NmrScreen', 'pH', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:16:06_00002'), 
 ('NMRS', 'NmrScreen', 'temperature', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:16:06_00001'), 
 ('NMRS', 'NmrScreen', 'userProtocolCode', 'www.ccpn.ac.uk_Fogh_2012-07-23-11:52:27_00002'), 
 ('NMRS', 'TrialSet', 'evaluateOnlyUnambiguous', 'www.ccpn.ac.uk_Fogh_2012-07-16-12:06:02_00001'), 
 ('NMRS', 'TrialSet', 'evaluateSingleResonance', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00008'), 
 ('NMRS', 'TrialSet', 'experimentWeights', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00006'), 
 ('NMRS', 'TrialSet', 'identifyAllosteric', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00007'), 
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
 ('ANAY', 'WindowPanel', 'useMultiplPeakLists', 'useMultiplePeakLists', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00013'), 
 ('CALC', 'ParameterGroup', 'datas', 'data', 'www.ccpn.ac.uk_Fogh_2011-10-11-16:36:23_00002'), 
 ('NMRS', 'Mixture', 'trial', 'trials', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00015'), 
 ('NMRS', 'RegionWeight', 'maxppm', 'maxPpm', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00005'), 
]

# ValueType changes
# change types are : 'ignore': do nothing, 'delay': available for calculation
# (prefix, typeName, elemName, action, newGuid, elemMap, valueTypeGuid)
typeChanges = [
]

# Different elements with matching qualifiedNames
# (element.qName, differentTags, oldGuid, newGuid
nameMatches = [
]

# Differences for matching elements, 
# excluding those where only names and/or valueTypes differ
# (oldElem.qName, newElem.name, oldGuid, newGuid, differentTags
allDiffs = [
 ('cambridge.WmsProtocol.InterfaceParameter.hicard', 'hicard', 'www.ccpn.ac.uk_Fogh_2011-03-22-17:23:24_00005', 'www.ccpn.ac.uk_Fogh_2011-03-22-17:23:24_00005', set(['defaultValue', 'documentation', 'locard'])), 
 ('cambridge.WmsProtocol.InterfaceParameter.locard', 'locard', 'www.ccpn.ac.uk_Fogh_2011-03-22-17:23:24_00006', 'www.ccpn.ac.uk_Fogh_2011-03-22-17:23:24_00006', set(['defaultValue', 'documentation'])), 
 ('ccp.nmr.NmrScreen.RegionWeight.weight', 'weight', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00006', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00006', set(['defaultValue'])), 
 ('ccpnmr.Analysis.PeakDrawMethod', 'PeakDrawMethod', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00002', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00002', set(['enumeration'])), 
 ('ccpnmr.Analysis.PeakFindVolumeMethod', 'PeakFindVolumeMethod', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00003', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00003', set(['enumeration'])), 
 ('ccpnmr.Analysis.SymbolStyle', 'SymbolStyle', 'www.ccpn.ac.uk_Fogh_2008-05-05-18:37:53_00002', 'www.ccpn.ac.uk_Fogh_2008-05-05-18:37:53_00002', set(['enumeration'])), 
]
