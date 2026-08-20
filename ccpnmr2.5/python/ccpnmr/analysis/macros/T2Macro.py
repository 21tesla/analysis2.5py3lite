from ccpnmr.analysis.core.DataAnalysisBasic import calcT2List

def calculateT2List(argServer):
  
  t1List = argServer.getMeasurementList(className='T1List')
  
  if not t1List: # None avail or aborted
    msg = 'Stopping: Script cancelled or no T1 List available.'
    argServer.showWarning(msg )
    return

  t1rhoList = argServer.getMeasurementList(className='T1RhoList')
  
  if not t1rhoList:
    msg = 'Stopping: Script cancelled or no T1rho List available.'
    argServer.showWarning(msg )
    return 
  
  spinLock = argServer.askFloat('Enter spin lock:', 0.0)
  
  if spinLock is None:
    msg = 'Stopping: Script cancelled.'
    argServer.showWarning(msg )
    return 
  
  # TBC the two lists could use different units
  # but highly unlikely to use anything other than seconds
  
  t2List = calcT2List(t1List, t1rhoList, spinLock)
  
  if t2List:
    argServer.parent.editMeasurements(t2List)
  
  return t2List

