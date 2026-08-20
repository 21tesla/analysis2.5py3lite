#
# Exercise 3.1: Create constraint list...
#
# - Distance (DO THIS AS EXAMPLE!)
# - Dihedral constraints
#

from memops.api import Implementation

if __name__ == '__main__':

  ccpnProject = Implementation.MemopsRoot(name = 'example3')
  nmrProject = ccpnProject.newNmrProject(name = 'test')
  
  nmrConstraintStore = ccpnProject.newNmrConstraintStore(nmrProject = nmrProject)
    
  #
  # Make some 'fixed' resonances - note that again these don't have to be linked
  #
  # The FixedAtomSet/FixedResonanceSet structure works exactly the same as in the
  # Nmr package, except that here it is in principle 'frozen' in time.
  #
  
  for i in range(10):
    nmrConstraintStore.newFixedResonance(isotopeCode = '1H')  
    nmrConstraintStore.newFixedResonance(isotopeCode = '13C')
    
  #
  # Nmr Distance Constraints
  #
  # Note that can set links to Nmr.Experiments and MeasurementLists that were relevant
  # for generating this constraint list!
  #
  
  dcl = nmrConstraintStore.newDistanceConstraintList(name = 'test', details = 'Fake list')
  
  #
  # Make constraints with constraint items. Note that in the case of distance constraints
  # these are a subclass of 'PairwiseConstraintItems', and the link to the resonances can
  # be 'ordered' for output purposes
  #
  # None of the below keywords are obligatory when setting (some are however set with default
  # values)
  #
  
  kywdsBasic  = {'details': 'Fake', 'origData': 534523.12, 'weight': 1}
  kywdsTarget = {'targetValue': 3.5, 'error': 0.5}
  kywdsUpLo   = {'lowerLimit':  2.0, 'upperLimit': 5.0}
  
  protonResonanceList = list(nmrConstraintStore.findAllFixedResonances(isotopeCode = '1H'))
  carbonResonanceList = list(nmrConstraintStore.findAllFixedResonances(isotopeCode = '13C'))
  
  for i in range(5):
  
    #
    # Make a distance constraint.
    #
    # Note: can also set values for backCalcVolume.
    #
    
    kywds = kywdsBasic.copy()
    
    if i % 2:
      kywds.update(kywdsTarget)
    else:
      kywds.update(kywdsUpLo)
  
    distanceConstraint = dcl.newDistanceConstraint(**kywds)
    
    for j in range(2):
      
      #
      # Make a distance constraint item, set so that 1H resonances always come 'first'.
      # Note that the resonances link is stored as a set, so the input order does *not* matter
      #
      
      index = i * 2 + j

      resonance1 = protonResonanceList[index - 1]
      resonance2 = carbonResonanceList[index - 1]
      
      distanceConstraint.newDistanceConstraintItem(resonances = [resonance2,resonance1], firstResonance = resonance1)
  
  
  #
  # Get the info out again - note the use of 'orderedResonances'. If firstResonance is *not* set this will return
  # an empty list!
  # 
  
  for dc in dcl.sortedConstraints():
  
    infoDict = {}

    for attrName in ('lowerLimit','upperLimit','targetValue'):
      value = getattr(dc,attrName)
      if value == None:
        value = 'n/a'
      else:
        value = "%6.2f" % value
      infoDict[attrName] = value
  
    print "%-4d   %-10s %-10s %-10s" % (dc.serial,infoDict['lowerLimit'],infoDict['targetValue'],infoDict['upperLimit'])

    for dci in dc.sortedItems():

      resonances = dci.orderedResonances
      if not resonances:
        resonances = dci.sortedResonances()

      print "  Resonances %3d - %3d" % (resonances[0].serial,resonances[1].serial)
      
    print


    
    
  #
  # Nmr Dihedral Constraints
  #
  
  dhcl = nmrConstraintStore.newDihedralConstraintList(name = 'test', details = 'Fake list')
  
  print dhcl.serial  # Note that this is 2 because distance list has serial 1!
  
  #
  # For dihedral constraints, the resonances are linked on the constraint level.
  #
  # Note that the link is ordered by default, so the order the resonances go in is the
  # order they come out. The output is a list, and sortedResonances() method does not exist.
  #
  
  for i in range(5):
  
    resonanceList = []
    
    #
    # Get H-C-C-H resonances out (because of way that was set up, in real life has to be
    # more sophisticated!)
    #
    
    resonanceList.append(nmrConstraintStore.findFirstFixedResonance(serial = i + 1))  # proton
    resonanceList.append(nmrConstraintStore.findFirstFixedResonance(serial = i + 2))  # carbon
    resonanceList.append(nmrConstraintStore.findFirstFixedResonance(serial = i + 4))  # carbon
    resonanceList.append(nmrConstraintStore.findFirstFixedResonance(serial = i + 3))  # proton
      
    dihedralConstraint = dhcl.newDihedralConstraint(resonances = resonanceList)
    
    #
    # Below constrains the dihedral angle from 100 to 140 degrees and -60 to 20 degrees
    #
    
    dihedralConstraintItem = dihedralConstraint.newDihedralConstraintItem(lowerLimit = 100.0, upperLimit = 140.0, targetValue = 130.0)
    dihedralConstraintItem = dihedralConstraint.newDihedralConstraintItem(lowerLimit = 300.0, upperLimit = 20.0,  targetValue = 0.0)

  for dhdlConst in dhcl.sortedConstraints():
    print 'Constraint [%d]: resonances [%s] - [%s] - [%s] - [%s]' % (
      dhdlConst.serial,
      dhdlConst.resonances[0].serial, dhdlConst.resonances[1].serial,
      dhdlConst.resonances[2].serial, dhdlConst.resonances[3].serial)

    for dhdlConstItem in dhdlConst.sortedItems():
      print 'Item [%d]: [%.2f] - [%.2f]; [%.2f]' % (
        dhdlConstItem.serial, dhdlConstItem.lowerLimit,
        dhdlConstItem.upperLimit, dhdlConstItem.targetValue)

    print
