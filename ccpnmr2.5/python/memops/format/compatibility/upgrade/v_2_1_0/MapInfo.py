
# Packages, classElements and AbstractDataTypes skipped in new model
# (prefix, typeName, elemName, newGuid, elemType)
skipElements = [
 ('IMPL', 'MemopsRoot', 'currentNmrScreenStore', 'ccpn_automatic_memops.Implementation.MemopsRoot.currentNmrScreenStore', 'MetaRole'), 
 ('IMPL', 'MemopsRoot', 'nmrScreenStores', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:34_00002', 'MetaRole'), 
 ('MOLE', 'Molecule', 'compounds', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:51:56_00027', 'MetaRole'), 
 ('NMRS', 'Compound', None, 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00003', 'MetaClass'), 
 ('NMRS', 'MixtureComponent', 'compound', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00003', 'MetaRole'), 
 ('NMRS', 'NmrScreen', 'nmrScreenStore', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00043', 'MetaRole'), 
 ('NMRS', 'NmrScreen', 'trials', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00029', 'MetaRole'), 
 ('NMRS', 'NmrScreenStore', None, 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00002', 'MetaClass'), 
 ('NMRS', 'RefSpectrum', None, 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00004', 'MetaClass'), 
 ('NMRS', 'Trial', 'nmrScreen', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00028', 'MetaRole'), 
 ('NMRS', 'Trial', 'trialExperiments', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:51:56_00006', 'MetaRole'), 
 ('NMRS', 'TrialExperiment', 'trial', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:51:56_00005', 'MetaRole'), 
 ('NMRS', 'TrialHit', 'refSpectrums', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00047', 'MetaRole'), 
]

# classElements skipped in new model, but available for simple data transfer
# (prefix, typeName, elemName, newGuid, elemMap, valueTypeGuid)
delayElements = [
 ('ANAY', 'WindowPanelGroup', 'showSameSpectra', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:04:35_00004', {'tag': 'ANAY.WindowPanelGroup.showSameSpectra', 'type': 'attr', 'name': 'showSameSpectra'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00028'), 
 ('CALC', 'NmrCalcStore', 'methodStoreName', 'www.ccpn.ac.uk_Fogh_2010-05-10-13:46:58_00003', {'eType': 'cplx', 'tag': 'CALC.NmrCalcStore.methodStoreName', 'type': 'attr', 'name': 'methodStoreName'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00033'), 
 ('CALC', 'NmrCalcStore', 'softwareName', 'www.ccpn.ac.uk_Fogh_2010-05-10-13:46:58_00004', {'eType': 'cplx', 'tag': 'CALC.NmrCalcStore.softwareName', 'type': 'attr', 'name': 'softwareName'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00033'), 
 ('CALC', 'NmrCalcStore', 'softwareVersion', 'www.ccpn.ac.uk_Fogh_2010-05-10-13:46:58_00005', {'eType': 'cplx', 'tag': 'CALC.NmrCalcStore.softwareVersion', 'type': 'attr', 'name': 'softwareVersion'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00033'), 
 ('NMRS', 'Mixture', 'pH', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00019', {'tag': 'NMRS.Mixture.pH', 'type': 'attr', 'name': 'pH'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:53_00031'), 
 ('NMRS', 'Mixture', 'solvent', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00020', {'eType': 'cplx', 'tag': 'NMRS.Mixture.solvent', 'type': 'attr', 'name': 'solvent'}, 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00005'), 
 ('NMRS', 'Mixture', 'volume', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00018', {'tag': 'NMRS.Mixture.volume', 'type': 'attr', 'name': 'volume'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:54_00007'), 
 ('NMRS', 'MixtureComponent', 'concentration', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00025', {'tag': 'NMRS.MixtureComponent.concentration', 'type': 'attr', 'name': 'concentration'}, 'www.ccpn.ac.uk_Fogh_2006-08-16-14:22:54_00007'), 
]

# MetaConstraints added in new model
# (qualifiedName, guid)
newConstraints = [
 ('cambridge.WmsProtocol.InterfaceParameter.hicard.hicard_consistent_with_ProtocolParameter_hicard', 'www.ccpn.ac.uk_Fogh_2013-10-11-09:59:51_00001'), 
 ('cambridge.WmsProtocol.InterfaceParameter.locard.locard_consistent_with_ProtocolParameter_locard', 'www.ccpn.ac.uk_Fogh_2013-10-11-09:59:51_00002'), 
 ('ccp.nmr.NmrCalc.Run.masterRun.derived_runs_cannot_be_nested', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00003'), 
 ('ccp.nmr.NmrReference.ChemAtomNmrDistrib.consistent_number_valuesPerPoint', 'www.ccpn.ac.uk_Fogh_2012-04-13-14:02:18_00001'), 
 ('ccp.nmr.NmrReference.ChemAtomNmrDistrib.nd_distribution_is_normalised', 'www.ccpn.ac.uk_Fogh_2012-04-13-14:02:18_00002'), 
 ('ccp.nmr.NmrReference.ChemAtomNmrDistrib.numbers_are_positive', 'www.ccpn.ac.uk_Fogh_2012-04-13-14:02:18_00003'), 
 ('ccp.nmr.NmrReference.ChemAtomNmrDistrib.refAtoms.len_refatoms_eq_ndim', 'www.ccpn.ac.uk_Fogh_2012-04-13-13:40:44_00005'), 
 ('ccp.nmr.NmrScreen.ExperimentHit.normalisedChange.absvalue_le_1', 'www.ccpn.ac.uk_Fogh_2012-04-18-15:31:23_00002'), 
 ('ccp.nmr.NmrScreen.RegionWeight.intervals_do_not_overlap', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:06_00005'), 
 ('ccp.nmr.NmrScreen.RegionWeight.minPpm_lt_maxPpm', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:06_00004'), 
]

# Mandatory classElements added in new model
# New ClassElements with locard !=0, no default, not derived or Implementation
# (prefix, typeName, elemName, newGuid)
newMandatories = [
 ('NMRR', 'ChemAtomNmrDistrib', 'valuesPerPoint', 'www.ccpn.ac.uk_Fogh_2012-04-13-13:40:44_00003'), 
 ('NMRS', 'ExperimentWeight', 'expCode', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00001'), 
 ('NMRS', 'ExperimentWeight', 'trialSet', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00005'), 
 ('NMRS', 'MixtureComponent', 'componentName', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00004'), 
 ('NMRS', 'MixtureComponent', 'componentType', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00005'), 
 ('NMRS', 'NmrScreen', 'memopsRoot', 'ccpn_automatic_ccp.nmr.NmrScreen.NmrScreen.memopsRoot'), 
 ('NMRS', 'RegionWeight', 'maxPpm', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00005'), 
 ('NMRS', 'RegionWeight', 'minPpm', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00004'), 
 ('NMRS', 'RegionWeight', 'serial', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00003'), 
 ('NMRS', 'RegionWeight', 'trialSet', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00010'), 
 ('NMRS', 'Trial', 'trialSet', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00008'), 
 ('NMRS', 'TrialExperiment', 'expCode', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00007'), 
 ('NMRS', 'TrialExperiment', 'mixture', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00001'), 
 ('NMRS', 'TrialHit', 'componentName', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00011'), 
 ('NMRS', 'TrialSet', 'nmrScreen', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00001'), 
 ('NMRS', 'TrialSet', 'serial', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00012'), 
 ('REFD', 'RefDataStore', 'memopsRoot', 'ccpn_automatic_ccp.lims.RefData.RefDataStore.memopsRoot'), 
 ('REFD', 'RefDataStore', 'refSampleComponentStore', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00019'), 
 ('REFD', 'RefNmrSpectrum', 'componentName', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00021'), 
 ('REFD', 'RefNmrSpectrum', 'refDataStore', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00016'), 
]

# Packages, classElements and AbstractDataTypes added in new model
# Optional, i.e. excluding mandatory classElements given above
# (prefix, typeName, elemName, newGuid)
newElements = [
 ('ANA3', 'AnalysisProjectV3', 'details', 'www.ccpn.ac.uk_Fogh_2012-04-13-13:40:40_00002'), 
 ('ANA3', 'SpectrumView', 'intensityScaling', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:11_00001'), 
 ('ANA3', 'SpectrumView', 'isInToolbar', 'www.ccpn.ac.uk_Fogh_2012-05-14-13:21:11_00001'), 
 ('ANAL', 'SpectrumWindow', 'isZeroLineShown', 'www.ccpn.ac.uk_Fogh_2013-05-03-11:50:01_00001'), 
 ('ANAL', 'SpectrumWindowView', 'isContourLineVisible', 'www.ccpn.ac.uk_Fogh_2013-05-07-17:07:06_00001'), 
 ('ANAW', 'AbstractModule', 'details', 'www.ccpn.ac.uk_Fogh_2012-09-10-14:34:37_00001'), 
 ('ANAW', 'AbstractModule', 'keywords', 'www.ccpn.ac.uk_Fogh_2012-09-10-14:34:37_00002'), 
 ('ANAY', 'AbstractPanel', 'layoutArea', 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:29_00001'), 
 ('ANAY', 'AbstractPanel', 'name', 'www.ccpn.ac.uk_Fogh_2012-05-03-14:05:48_00001'), 
 ('ANAY', 'AbstractPanel', 'rank', 'www.ccpn.ac.uk_Fogh_2012-04-18-15:31:21_00001'), 
 ('ANAY', 'Layout', 'keywords', 'www.ccpn.ac.uk_Fogh_2012-05-14-13:21:09_00001'), 
 ('ANAY', 'Layout', 'rank', 'www.ccpn.ac.uk_Fogh_2012-09-10-14:34:35_00001'), 
 ('ANAY', 'LayoutArea', None, 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:27_00002'), 
 ('ANAY', 'SpectrumSharing', None, 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:27_00001'), 
 ('ANAY', 'WindowPanel1d', 'showPeakPickLevel', 'www.ccpn.ac.uk_Fogh_2012-09-17-10:38:43_00001'), 
 ('ANAY', 'WindowPanelGroup', 'gridCell', 'www.ccpn.ac.uk_Fogh_2012-04-13-13:34:58_00002'), 
 ('ANAY', 'WindowPanelGroup', 'gridSpan', 'www.ccpn.ac.uk_Fogh_2012-04-13-13:40:40_00001'), 
 ('ANAY', 'WindowPanelGroup', 'isGridGroup', 'www.ccpn.ac.uk_Fogh_2012-04-13-13:34:58_00001'), 
 ('ANAY', 'WindowPanelGroup', 'spectrumSharing', 'www.ccpn.ac.uk_Fogh_2012-08-16-17:30:29_00002'), 
 ('CALC', 'Run', 'derivedRuns', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00001'), 
 ('CALC', 'Run', 'masterRun', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00002'), 
 ('CALC', 'Run', 'methodStoreName', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00004'), 
 ('CALC', 'Run', 'softwareName', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00005'), 
 ('CALC', 'Run', 'softwareVersion', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00006'), 
 ('CALC', 'Run', 'wmsProtocolName', 'www.ccpn.ac.uk_Fogh_2012-06-04-14:36:41_00007'), 
 ('IMPL', 'MemopsRoot', 'currentNmrScreen', 'ccpn_automatic_memops.Implementation.MemopsRoot.currentNmrScreen'), 
 ('IMPL', 'MemopsRoot', 'currentRefDataStore', 'ccpn_automatic_memops.Implementation.MemopsRoot.currentRefDataStore'), 
 ('IMPL', 'MemopsRoot', 'nmrScreens', 'ccpn_automatic_memops.Implementation.MemopsRoot.nmrScreen'), 
 ('IMPL', 'MemopsRoot', 'refDataStores', 'ccpn_automatic_memops.Implementation.MemopsRoot.refDataStore'), 
 ('NMR', 'Experiment', 'userExpCode', 'www.ccpn.ac.uk_Fogh_2012-07-23-11:52:27_00001'), 
 ('NMR', 'Peak', 'height', 'www.ccpn.ac.uk_Fogh_2012-04-12-17:58:12_00001'), 
 ('NMR', 'Peak', 'volume', 'www.ccpn.ac.uk_Fogh_2012-04-12-17:58:12_00002'), 
 ('NMRC', 'FixedResonance', 'covalentlyBound', 'www.ccpn.ac.uk_Fogh_2012-06-25-14:41:56_00001'), 
 ('NMRS', 'ComponentType', None, 'www.ccpn.ac.uk_Fogh_2012-03-28-17:18:50_00001'), 
 ('NMRS', 'ExpCode', None, 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:06_00001'), 
 ('NMRS', 'ExperimentHit', 'normalisedChange', 'www.ccpn.ac.uk_Fogh_2012-04-18-15:31:23_00001'), 
 ('NMRS', 'ExperimentWeight', None, 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:44_00001'), 
 ('NMRS', 'ExperimentWeight', 'changeIsPositive', 'www.ccpn.ac.uk_Fogh_2012-07-09-10:23:18_00001'), 
 ('NMRS', 'ExperimentWeight', 'intensityScale', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00003'), 
 ('NMRS', 'ExperimentWeight', 'meritThreshold', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00002'), 
 ('NMRS', 'ExperimentWeight', 'weight', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00004'), 
 ('NMRS', 'Mixture', 'trialExperiments', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00002'), 
 ('NMRS', 'NmrScreen', 'pH', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:16:06_00002'), 
 ('NMRS', 'NmrScreen', 'refDataStoreNames', 'www.ccpn.ac.uk_Fogh_2012-04-12-17:58:11_00001'), 
 ('NMRS', 'NmrScreen', 'temperature', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:16:06_00001'), 
 ('NMRS', 'NmrScreen', 'trialSets', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00002'), 
 ('NMRS', 'NmrScreen', 'userProtocolCode', 'www.ccpn.ac.uk_Fogh_2012-07-23-11:52:27_00002'), 
 ('NMRS', 'RegionWeight', None, 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:06_00003'), 
 ('NMRS', 'RegionWeight', 'weight', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00006'), 
 ('NMRS', 'TrialHit', 'isConfirmed', 'www.ccpn.ac.uk_Fogh_2012-05-18-12:04:57_00001'), 
 ('NMRS', 'TrialHit', 'refNmrSpectra', 'www.ccpn.ac.uk_Fogh_2012-03-29-15:58:32_00001'), 
 ('NMRS', 'TrialSet', None, 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:06_00002'), 
 ('NMRS', 'TrialSet', 'details', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00016'), 
 ('NMRS', 'TrialSet', 'evaluateOnlyUnambiguous', 'www.ccpn.ac.uk_Fogh_2012-07-16-12:06:02_00001'), 
 ('NMRS', 'TrialSet', 'evaluateSingleResonance', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00008'), 
 ('NMRS', 'TrialSet', 'experimentWeights', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00006'), 
 ('NMRS', 'TrialSet', 'identifyAllosteric', 'www.ccpn.ac.uk_Fogh_2012-07-06-13:03:50_00007'), 
 ('NMRS', 'TrialSet', 'name', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00013'), 
 ('NMRS', 'TrialSet', 'regionWeights', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00011'), 
 ('NMRS', 'TrialSet', 'trials', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00009'), 
 ('NMRS', 'TrialSet', 'useInverseEffects', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00014'), 
 ('NMRS', 'TrialSet', 'useVolume', 'www.ccpn.ac.uk_Fogh_2012-05-21-18:09:12_00015'), 
 ('REFD', None, None, 'www.ccpn.ac.uk_Fogh_2012-03-28-17:19:49_00001'), 
 ('REFD', 'RefDataStore', None, 'www.ccpn.ac.uk_Fogh_2012-03-28-17:19:49_00002'), 
 ('REFD', 'RefDataStore', 'name', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00020'), 
 ('REFD', 'RefDataStore', 'refNmrSpectra', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00017'), 
 ('REFD', 'RefNmrSpectrum', None, 'www.ccpn.ac.uk_Fogh_2012-03-29-14:18:52_00001'), 
 ('REFD', 'RefNmrSpectrum', 'temperature', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00022'), 
 ('REFD', 'RefNmrSpectrum', 'trialHits', 'www.ccpn.ac.uk_Fogh_2012-03-29-15:58:32_00002'), 
 ('REFS', 'RefSampleComponentStore', 'refDataStores', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00018'), 
 ('SAM', 'AbstractSample', 'solvent', 'www.ccpn.ac.uk_Fogh_2012-03-28-17:22:44_00014'), 
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
 ('ANAY', 'WindowPanel', 'panelAxiss', 'panelAxes', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00009'), 
 ('ANAY', 'WindowPanel', 'useMultiplPeakLists', 'useMultiplePeakLists', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00013'), 
 ('ANAY', 'WindowPanel', 'windowPanelGroups', 'windowPanelGroup', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:04:12_00002'), 
 ('ANAY', 'WindowPanel1D', None, 'WindowPanel1d', 'www.ccpn.ac.uk_Fogh_2011-11-16-17:07:15_00023'), 
 ('ANAY', 'WindowPanel1D', 'labelAngle', 'WindowPanel1d', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:04:12_00001'), 
 ('ANAY', 'WindowPanel1D', 'showAsStack', 'WindowPanel1d', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00021'), 
 ('ANAY', 'WindowPanel1D', 'showIntegrals', 'WindowPanel1d', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00020'), 
 ('ANAY', 'WindowPanel1D', 'stackOffset', 'WindowPanel1d', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00022'), 
 ('ANAY', 'WindowPanel1D', 'useAutoScale', 'WindowPanel1d', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00019'), 
 ('ANAY', 'WindowPanel1D', 'valueScale', 'WindowPanel1d', 'www.ccpn.ac.uk_Fogh_2011-11-30-11:03:44_00018'), 
 ('CALC', 'ParameterGroup', 'datas', 'data', 'www.ccpn.ac.uk_Fogh_2011-10-11-16:36:23_00002'), 
 ('NMRS', 'Trial', 'mixture', 'mixtures', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00014'), 
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
 ('ccp.lims.Holder.Holder.holderCategories', 'holderCategories', 'www.ccpn.ac.uk_Fogh_2006-08-16-18:23:27_00006', 'www.ccpn.ac.uk_Fogh_2006-08-16-18:23:27_00006', set(['locard'])), 
 ('ccp.lims.Sample.AbstractSample.sampleCategories', 'sampleCategories', 'www.ccpn.ac.uk_Fogh_2006-08-16-18:22:46_00003', 'www.ccpn.ac.uk_Fogh_2006-08-16-18:22:46_00003', set(['locard'])), 
 ('ccp.nmr.NmrReference.ChemAtomNmrDistrib.refAtoms', 'refAtoms', 'www.ccpn.ac.uk_Fogh_2010-05-14-17:17:46_00001', 'www.ccpn.ac.uk_Fogh_2010-05-14-17:17:46_00001', set(['changeability'])), 
 ('ccp.nmr.NmrScreen.NmrScreen', 'NmrScreen', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00006', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00006', set(['supertype', 'parentRole', 'partitionsChildren', 'supertypes'])), 
 ('ccp.nmr.NmrScreen.NmrScreen.startDate', 'startDate', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00038', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00038', set(['locard'])), 
 ('ccp.nmr.NmrScreen.RefSpectrum.concentration', 'concentration', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00053', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00053', set(['container'])), 
 ('ccp.nmr.NmrScreen.RefSpectrum.dataSource', 'dataSource', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00050', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00050', set(['container'])), 
 ('ccp.nmr.NmrScreen.RefSpectrum.details', 'details', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00056', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00056', set(['container'])), 
 ('ccp.nmr.NmrScreen.RefSpectrum.name', 'name', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00052', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00052', set(['container'])), 
 ('ccp.nmr.NmrScreen.RefSpectrum.pH', 'pH', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00054', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00054', set(['container'])), 
 ('ccp.nmr.NmrScreen.RefSpectrum.serial', 'serial', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00051', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00051', set(['container'])), 
 ('ccp.nmr.NmrScreen.RefSpectrum.solvent', 'solvent', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00055', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00055', set(['container'])), 
 ('ccp.nmr.NmrScreen.Solvent', 'Solvent', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00005', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00005', set(['container', 'enumeration'])), 
 ('ccp.nmr.NmrScreen.Trial', 'Trial', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00009', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00009', set(['parentRole'])), 
 ('ccp.nmr.NmrScreen.Trial.mixture', 'mixtures', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00014', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:38_00014', set(['hicard', 'name', 'locard'])), 
 ('ccp.nmr.NmrScreen.TrialExperiment', 'TrialExperiment', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00012', 'www.ccpn.ac.uk_Fogh_2009-11-19-14:50:32_00012', set(['parentRole'])), 
 ('ccpnmr.Analysis.PeakDrawMethod', 'PeakDrawMethod', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00002', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00002', set(['enumeration'])), 
 ('ccpnmr.Analysis.PeakFindVolumeMethod', 'PeakFindVolumeMethod', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00003', 'www.ccpn.ac.uk_Fogh_2006-10-03-11:26:03_00003', set(['enumeration'])), 
 ('ccpnmr.Analysis.SymbolStyle', 'SymbolStyle', 'www.ccpn.ac.uk_Fogh_2008-05-05-18:37:53_00002', 'www.ccpn.ac.uk_Fogh_2008-05-05-18:37:53_00002', set(['enumeration'])), 
]
