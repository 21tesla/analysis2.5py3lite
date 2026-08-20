#
# Exercise 2.1: Create some measurements
#
# - Shifts (DO THIS AS EXAMPLE!)
# - JCoupling
#

from memops.api import Implementation

if __name__ == '__main__':

  ccpnProject = Implementation.MemopsRoot(name = 'example2')
  nmrProject = ccpnProject.newNmrProject(name = 'test')
    
  #
  # Make some resonances - note that these don't have to be linked
  #
  
  for i in range(10):
    nmrProject.newResonance(isotopeCode = '1H')  

  #
  # Nmr measurements - Shifts
  #
    
  #
  # Create a shift list
  #
  
  shiftList = nmrProject.newShiftList(name = 'My new list', details = 'Fake information')
  
  #
  # Make shifts... note that a Resonance can only have 1 Shift value for a particular ShiftList!
  #
  # For Shift creation, one can also do:
  #
  # Nmr.Shift(shiftList, ...)
  #
  
  for i in range(10):
    resonance = nmrProject.findFirstResonance(serial = i+1)
    shiftList.newShift(resonance = resonance, value = i, error = 0.01, figOfMerit = 0.0, details = 'Fake')
    
  #
  # Get info back out, via Resonances (can also do this via Shifts of course)
  #
  
  print "Shiftlist information:"
  print
  
  shiftUnit = shiftList.unit
  
  for resonance in nmrProject.sortedResonances():   # If use .resonances, order is random
    shift = resonance.findFirstShift(parent = shiftList)
    
    if shift:
      print "Resonance %3d: %7.3f %s +/- %4.3f %s     %s" % (resonance.serial,shift.value,shiftUnit,shift.error,shiftUnit,shift.details)
  
  print
  print
  
  
  #
  # Nmr measurements - J couplings
  #
  
  jCouplingList = nmrProject.newJCouplingList()
  
  for i in range(5):
    resonance1 = nmrProject.findFirstResonance(serial = i+1)
    resonance2 = nmrProject.findFirstResonance(serial = i+6)
    
    jCoupling = jCouplingList.newJCoupling(resonances = [resonance1,resonance2], value = i)
    
  #
  # Get info back out...
  #
  # Note the use of measurement - this is because the NMR measurement lists and measurement classes
  # are all subclasses of abstract supersets. This is sometimes (but not always?) reflected in the
  # link names.
  #
  
  print "J coupling list information:"
  print
  
  for jCoupling in jCouplingList.sortedMeasurements():
  
    resonances = jCoupling.sortedResonances()
    print "Resonances %3d - %3d:  %7.3f" % (resonances[0].serial,resonances[1].serial,jCoupling.value)
  
  
