
from ccp.util import NmrCalc as nmrCalcUtil
from ccp.api.nmr import NmrCalc

ISD = 'ISD'
STRING = 'STRING'
INT = 'INT'
FLOAT = 'FLOAT'
BOOL = 'BOOL'
INTLIST = 'INTLIST'
STRINGLIST = 'INTLIST'
STRINGDICT = 'STRINGDICT'

PROVISIONAL = 'provisional'

translate32charLimit = {'ccpn.export.ms_name':'ccpn.export.molecular_system_name'}
 
SIM_ATTRS = [('name',STRING),
             ('working_path',STRING),
             ('temp_path',STRING),
             ('shared_temp_path',BOOL),
             ('cns_executable',STRING),
             ('fileroot',STRING),
             ('temperature',FLOAT),
             ('naming_system',STRING),
             ('pdb_files.ensemble',STRING),
             ('pdb_files.n',INT),
             ('pdb_files.residue_list',INTLIST),
             ('ccpn.project_filename',STRING),
             ('ccpn.nmr_project_name',STRING),
             ('ccpn.export.ms_name',STRING),
             ('ccpn.export.project.key',STRING),
             ('ccpn.export.project.enabled',BOOL),
             ('ccpn.export.ensemble.enabled',BOOL),
             ('molecule.filename',STRING),
             ('molecule.format',STRING),
             ('molecule.key',STRING),
             ('molecule.first_residue_number',INT),
             ('molecule.initial_conformation',STRING),
             ('molecule.exclude_hydrogens',BOOL),
             ('replica.n_replicas',INT),
             ('replica.n_samples',INT),
             ('replica.hmc.steps',INT),
             ('replica.hmc.md.steps',INT),
             ('replica.hmc.stepsize',FLOAT),
             ('replica.hmc.adjust_stepsize',FLOAT),
             ('replica.communication',STRING),
             ('replica.override_parameters',BOOL),
             ('replica.save_interval',INT),
             ('replica.full_save',BOOL),
             ('replica.python',STRING),
             ('replica.niceness',INT),
             ('replica.background',BOOL),
             ('replica.name_server',STRING),
             ('replica.host_list',STRINGLIST),
             ('replica.temp_paths',STRINGDICT),
             ('replica.weight_schedule.initial',FLOAT),
             ('replica.weight_schedule.final',FLOAT),
             ('replica.weight_schedule.slope',FLOAT),
             ('replica.weight_schedule.first',INT),
             ('replica.weight_schedule.last',INT),
             ('replica.prior_schedule.initial',FLOAT),
             ('replica.prior_schedule.final',FLOAT),
             ('replica.prior_schedule.slope',FLOAT),
             ('replica.prior_schedule.first',INT),
             ('replica.prior_schedule.last',INT),
             ('analysis.burnin',INT),
             ('analysis.report.auto',BOOL),
             ('analysis.report.keep_sources',BOOL),
             ('analysis.pdb_viewer',STRING),
             ('analysis.gnuplot',STRING),
             ('analysis.pdf_latex',STRING),
             ('analysis.eps_to_pdf',STRING),
             ('analysis.eps_to_eps',STRING),
             ('analysis.max_samples',INT),
             ('analysis.whatif.binary',STRING),
             ('analysis.whatif.use_whatcheck',BOOL),
             ('analysis.whatif.show_traces',BOOL),
             ('analysis.whatif.enabled',BOOL),
             ('analysis.procheck.binary',STRING),
             ('analysis.procheck.show_traces',BOOL),
             ('analysis.procheck.enabled',BOOL),
             ('analysis.dssp.binary',STRING),
             ('analysis.dssp.enabled',BOOL),
             ('use_xterm',BOOL)]
            
EXP_DATA_GEN_ATTRS = [('data_type',STRING),
                      ('filename',STRING),
                      ('format',STRING),
                      ('key',STRING),
                      ('name',STRING),]

EXP_DATA_ATTR_DICT = {'noesy':(('theory.scale.update',BOOL),
                            ('error_model.error.initial',FLOAT),
                            ('error_model.error.update',BOOL)),
                   'distance':(('theory.scale.update',BOOL),
                            ('error_model.error.initial',FLOAT),
                            ('error_model.error.update',BOOL)),
                   'jcoupling':(('theory.karplus_curve.A',FLOAT),
                                ('theory.karplus_curve.B',FLOAT),
                                ('theory.karplus_curve.C',FLOAT),
                                ('theory.karplus_curve.update',BOOL),
                                ('error_model.error.initial',FLOAT),
                                ('error_model.error.update',BOOL)),
                   'rdc':(('theory.saupe_tensor.s1',FLOAT),
                          ('theory.saupe_tensor.s2',FLOAT),
                          ('theory.saupe_tensor.s3',FLOAT),
                          ('theory.saupe_tensor.s4',FLOAT),
                          ('theory.saupe_tensor.s5',FLOAT),
                          ('theory.saupe_tensor.update',BOOL),
                          ('error_model.error.initial',FLOAT),
                          ('error_model.error.update',BOOL)),
                   'dihedral':(('error_model.error.update',BOOL),),
                   'hbond':(('error_model.error.initial',FLOAT),
                            ('error_model.error.update',BOOL)),
                   'disulfide':(('distance',FLOAT),
                                ('error_model.error.update',BOOL),
                                ('error_model.error.initial',FLOAT))}
                                  
GET_RUN_PARAM_FUNCS = {INT:nmrCalcUtil.getRunIntParameter,
                       FLOAT:nmrCalcUtil.getRunFloatParameter,
                       BOOL:nmrCalcUtil.getRunBooleanParameter,
                       STRING:nmrCalcUtil.getRunTextParameter,
                       INTLIST:nmrCalcUtil.getRunIntParameterList,
                       STRINGLIST:nmrCalcUtil.getRunTextParameterList,
                       STRINGDICT:nmrCalcUtil.getRunTextParameterDict}
 
SET_RUN_PARAM_FUNCS = {INT:nmrCalcUtil.setRunIntParameter,
                       FLOAT:nmrCalcUtil.setRunFloatParameter,
                       BOOL:nmrCalcUtil.setRunBooleanParameter,
                       STRING:nmrCalcUtil.setRunTextParameter,
                       INTLIST:nmrCalcUtil.setRunIntParameterList,
                       STRINGLIST:nmrCalcUtil.setRunTextParameterList,
                       STRINGDICT:nmrCalcUtil.setRunTextParameterDict}

SET_DATA_PARAM_FUNCS = {INT:nmrCalcUtil.setDataIntParameter,
                        FLOAT:nmrCalcUtil.setDataFloatParameter,
                        BOOL:nmrCalcUtil.setDataBooleanParameter,
                        STRING:nmrCalcUtil.setDataTextParameter,
                        INTLIST:nmrCalcUtil.setDataIntParameterList,
                        STRINGLIST:nmrCalcUtil.setDataTextParameterList,
                        STRINGDICT:nmrCalcUtil.setDataTextParameterDict}

GET_DATA_PARAM_FUNCS = {INT:nmrCalcUtil.getDataIntParameter,
                        FLOAT:nmrCalcUtil.getDataFloatParameter,
                        BOOL:nmrCalcUtil.getDataBooleanParameter,
                        STRING:nmrCalcUtil.getDataTextParameter,
                        INTLIST:nmrCalcUtil.getDataIntParameterList,
                        STRINGLIST:nmrCalcUtil.getDataTextParameterList,
                        STRINGDICT:nmrCalcUtil.getDataTextParameterDict}
                                
def nmrCalcRunToIsd(run):
  from cambridge.isd.isd_project_template import sim
  from Isd.setup import setup_data, CCPN
  from CCPNReader import getObjectKeyString

  sim.data_sets = []
  sim.isd_version = 1.2 # Will be filled in properly after data model change
  sim.ccpn.nmr_project_name = run.nmrCalcStore.nmrProject.name
  sim.ccpn.export.project.key = getObjectKeyString(run)
                  
  for attrName, attrType in SIM_ATTRS:
  
    value = GET_RUN_PARAM_FUNCS[attrType](run, attrName)
    attrName = translate32charLimit.get(attrName, attrName)
    
    if value is not None:
      dotPath = attrName.split('.')
      
      subAttrName = dotPath.pop(0)
      obj = sim
      
      while dotPath:
        obj = getattr(obj, subAttrName)
        subAttrName = dotPath.pop(0)
    
      setattr(obj, subAttrName, value)
  
  
  for datum in run.inputs:
    if isinstance(datum, NmrCalc.ExternalData):
      nmrObjs = [None,]
      
    elif isinstance(datum, NmrCalc.ConstraintStoreData):
      nmrObjs = datum.constraintLists # Extant lists

    elif isinstance(datum, NmrCalc.MolSystemData):
      chains = datum.chains
      
      if chains:
        msCode = datum.molSystemCode
        chainCode = chains[0].code
        sim.ccpn.export.molecular_system_name = msCode
        sim.molecule.key = '%s|%s' % (msCode,chainCode)
        sim.molecule.format = CCPN
      
      continue
  
    else:
      continue
  
    values = []
    for attrName, attrType in EXP_DATA_GEN_ATTRS:
      value = GET_DATA_PARAM_FUNCS[attrType](datum, attrName)
      values.append(value)
    
    dataType, dataFileName, dataFormat, dataKey, dataName = values
    
    for nmrObj in nmrObjs:
      
      if nmrObj:
        key = getObjectKeyString(nmrObj)
        format = CCPN
        fileName = ''
      else:
        key = dataKey
        format = dataFormat
        fileName = dataFileName

      specs = setup_data(dataType, fileName, dataName, key, format)
      
      for attrName, attrType in EXP_DATA_ATTR_DICT[dataType]:
        value = GET_DATA_PARAM_FUNCS[attrType](datum, attrName)
 
        if value is not None:
          dotPath = attrName.split('.')
 
          subAttrName = dotPath.pop(0)
          obj = specs
 
          while dotPath:
            obj = getattr(obj, subAttrName)
            subAttrName = dotPath.pop(0)
 
          setattr(obj, subAttrName, value)
          
      sim.data_sets.append(specs)
 
  return sim
  
def getIsdNmrCalcStore(ccpnProject):

  nmrProject = ccpnProject.currentNmrProject
  calcStore = ccpnProject.findFirstNmrCalcStore(name=ISD, nmrProject=nmrProject) or \
              ccpnProject.newNmrCalcStore(name=ISD, nmrProject=nmrProject)
  
  return calcStore

 
def isdToNmrCalcRun(sim, ccpnProject, run=None):
 
  from Isd.setup import DISULFIDE, CCPN
  from CCPNReader import getKeysFromString, getObjectKeyString
 
  if run:
    calcStore = run.nmrCalcStore
 
    for datum in run.inputs:
      for param in datum.runParameters:
        param.delete()
      datum.delete()
 
  else:
    calcStore = getIsdNmrCalcStore(ccpnProject)
    runs = [r for r in calcStore.sortedRuns() if r.status == PROVISIONAL]
 
    if runs:
      run = runs[-1]
      
      for datum in run.input:
        for param in datum.runParameters:
          param.delete()
        datum.delete()

    else:
      run = calcStore.newRun(status=PROVISIONAL)
 
  nmrProject = calcStore.nmrProject
  sim.ccpn.nmr_project_name = nmrProject.name
  sim.ccpn.export.project.key = getObjectKeyString(run)
 
  for attrName1, attrType in SIM_ATTRS:
    attrName = translate32charLimit.get(attrName1, attrName1)
 
    dotPath = attrName.split('.')
 
    subAttrName = dotPath.pop(0)
    value = getattr(sim, subAttrName)
 
    while dotPath:
      subAttrName = dotPath.pop(0)
      value = getattr(value, subAttrName)
 
    SET_RUN_PARAM_FUNCS[attrType](run, attrName1, value)
 
  if (sim.molecule.format) == CCPN and sim.molecule.key:
    msCode, chainCode = sim.molecule.key.split('|')
    molSystem = ccpnProject.findFirstMolSystem(code=msCode)
 
    if molSystem:
      chain = molSystem.findFirstChain(code=chainCode)
 
      if chain:
        run.newMolSystemData(molSystemCode=msCode, chainCodes=[chainCode,])
 
  getCStore = ccpnProject.findFirstNmrConstraintStore
 
  for specs in sim.data_sets:
     dataType = specs.data_type
     dataName = specs.name
 
     if specs.format == CCPN:
       if dataType == DISULFIDE: # Not supported, yet
         continue
 
       csSerial, clSerial = getKeysFromString(specs.key)

       cStore = getCStore(serial=csSerial, nmrProject=nmrProject)
       if not cStore:
         continue
 
       constraintList = cStore.findFirstConstraintList(serial=clSerial)
       if not constraintList:
         continue

       if not dataName:
         dataName = '%s %d:%d' % (constraintList.className[:-14],
                                  csSerial, clSerial)
 
       datum = run.newConstraintStoreData(constraintStoreSerial=csSerial,
                                          constraintListSerials=[clSerial,],
                                          name=dataName)

     else:
       if not dataName:
         dataName = '%s format' % specs.format
 
       datum = run.newExternalData(name=dataName)

     for attrName, attrType in EXP_DATA_GEN_ATTRS:
       value = getattr(specs, attrName) or ''
       SET_DATA_PARAM_FUNCS[attrType](datum, attrName, value)
 
     for attrName, attrType in EXP_DATA_ATTR_DICT[dataType]:       
       dotPath = attrName.split('.')
       subAttrName = dotPath.pop(0)
       value = getattr(specs, subAttrName)
 
       while dotPath:
         subAttrName = dotPath.pop(0)
         value = getattr(value, subAttrName)
 
       SET_DATA_PARAM_FUNCS[attrType](datum, attrName, value)
 
  return run
