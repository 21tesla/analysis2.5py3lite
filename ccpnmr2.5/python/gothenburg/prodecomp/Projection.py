""" CCPN code for handling projection spectra.

NB this file does *not* import or use numpy, so the code can be freely used
"""
from memops.general.Implementation import ApiError

def getProjectionSpectra(nmrProject):
  
  spectra = []

  for experiment in nmrProject.sortedExperiments():
    acquExpDim = experiment.findFirstExpDim(isAcquisition=True)
    if not acquExpDim:
      experiment.sortedExpDims()[0].isAcquisition=True
    
    for dataSource in experiment.sortedDataSources():
      
      # filter for appropriate projection spectra  
      if dataSource.numDim != 2:
        # 2D exps only
        continue
      
      if dataSource.dataType != 'processed':
        #processed exps only
        continue
      
      dataDims = dataSource.sortedDataDims()
      dd0 = dataDims[0]
      shiftRefs = dd0.expDim.findAllExpDimRefs(measurementType='Shift')
      if len(shiftRefs) != 1:
        # only one shift reference in acq dimension
        continue

      if dd0.findFirstDataDimRef(expDimRef=list(shiftRefs)[0]) is None:
        continue
      
      if not dd0.expDim.isAcquisition:
        # First dim must be acquisition
        continue
      
      expDim1 = dataDims[1].expDim
      
      expDimRefs = expDim1.findAllExpDimRefs(measurementType='Shift')
      if len(expDimRefs) < 2:
        # second dim must have several Shift expDimRefs
        continue
      
      # displayNames for relevant ExpDimRefs must be unique
      names = [x.displayName for x in expDimRefs]
      if len(names) != len(set(names)):
        continue
      
      if projectionDimScalings(dataDims[1]):
        spectra.append(dataSource)

  return spectra

def projectionDimScalings(dataDim):
  """ Get DimensionScalings appropriate for projection spectra
  """
  result = []
  
  for ds in dataDim.sortedDimensionScalings():
    expDimRef =  ds.expDimRef
    if (len(ds.scalingFactors) == 1 and expDimRef.displayName is not None and
        expDimRef.measurementType=='Shift'):
      result.append(ds)
  #
  return result


def getIndirectShapeNames(dataSources):
  """Get names of all indirect shapes in Experiment and ExpDimref order
  """
  result = []
  
  experiments = []
  for dataSource in dataSources:
    experiment = dataSource.experiment
    if experiment not in experiments:
      experiments.append(experiment)
  
  for experiment in experiments:
    for expDim in experiment.sortedExpDims()[1:]:
      for expDimRef in expDim.sortedExpDimRefs():
        name = expDimRef.displayName
        if (name is not None and expDimRef.measurementType=='Shift'
            and name not in result):
          result.append(name)
  #
  return result
  
def getProjectionData(dataSources):
  """get list of shape names excluding acquisition dimension
  and spectrum,shapename scalingFactor matrix
  """
  
  allShapeNames = getIndirectShapeNames(dataSources)
  
  shapeNameSet = set()
  scalingFacs = []
  
  for dataSource in dataSources:
    dd = {}
    scalingFacs.append(dd)
    
    for dataDim in dataSource.sortedDataDims()[1:]:
      for dsc in projectionDimScalings(dataDim):
 
        sname = dsc.expDimRef.displayName
        if sname in dd:
          raise ApiError("%s ExpDimRef has duplicate displayName" % dsc)
 
        shapeNameSet.add(sname)
        dd[sname] = dsc.scalingFactors[0]
      
  # create defs matrix. NBNB currently must be int - may change later?
  # first get shapeNames in actual use in predetermined order
  shapeNames = [x for x in allShapeNames if x in shapeNameSet]
  defsMatrix = []
  nShapes = len(shapeNames)
  for dd in scalingFacs:
    ll = [int(dd.get(sname,0.0)) for sname in shapeNames]
    defsMatrix.append(ll)
  #
  return shapeNames, defsMatrix

def formatScalingFactor(factor):
  """ do string format of scaling factor
  """
  if factor == 0.0:
    x = None
  
  elif int(factor) == factor:
    x = '%+d' % int(factor)
      
  else:
    x = '%.3f' % factor
  #
  return x
