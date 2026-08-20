# MemopsRoot

key_project_settings = 'project.settings'

def get_isd_name():
  return 'ISD'

def print_distance_restraint(r):

  print '  number=%d, distance=%s, upper=%.2f, lower=%.2f' % \
        (r.number, str(r.distance), r.upper, r.lower)

  for a, b in r.contributions:
    print ' ',a, b

def getObjectKeyString(object, delimiter='|'):
  """Descrn: Make an object identifier string for a CCPN data model object
     Inputs: CCPN data model object
     Output: String
  """

  keys = object.getExpandedKey()
  
  for i in range(len(keys)):
    key = keys[i]
    
    keyType = type(key)
    if keyType is type([]):
      keys[i] = delimiter.join([str(k) for k in key])
    elif keyType is not type(''):
      keys[i] = str(key)
  
  return delimiter.join(keys)

def is_number(s):

  from string import digits

  for x in s:
    if not x in digits:
      return False

  return True

def getKeysFromString(word, delimiter='|'):
  """Descrn: Get a list of CCPN data model keys given an object identifier string. 
     Inputs: String
     Output: List of Keys (Words or Ints)
  """

  items = word.split(delimiter)
   
  keys = []
  for item in items:
    if is_number(item):
      key = int(item)
      
    else:
      key = item
      
    keys.append(key)  

  return keys

def get_ccpn_chain(ccp_project, key):

  """Descrn: Fetch a CCPN chain from a CCPN project using a list of object keys
     Inputs: Implementation.Project, List of keys (Words or Ints)
     Output: ccp.molecule.MolSystem.Chan
  """

  molSystemCode, chainCode = key
  
  molSystem = ccp_project.findFirstMolSystem(code=molSystemCode)
  
  if not molSystem:
    ValueError('No molecular system with code "%s" in CCPN project' % molSystemCode)
    
  chain = molSystem.findFirstChain(code=chainCode)
  
  if not chain:
    ValueError('No chain found with code "%s" in molecular system "%s"' % (chainCode, molSystemCode))

  return chain
    
def make_isd_atom(ccp_atom):
  """Descrn: Make an ISD Atom given a CCPN Atom object
             Also sets up the hetero atom name.
     Inputs: ccp.molecule.MolSystem.Atom
     Output: ISD Atom
  """

  from Atom import Atom
  
  ## add properties like mass, charge etc
    
  properties = {}

  isd_atom = Atom(ccp_atom.name, properties=properties)

  return isd_atom

def make_isd_residue(ccp_residue, connectivity):
  """Descrn: Make an ISD Residue given a CCPN Residue object
     Inputs: ccp.molecule.MolSystem.Residue
     Output: ISD residue
  """

  from Atom import Atom
  from ResidueIsd import Residue
        
  residue_type = ccp_residue.ccpCode
  residue_number = ccp_residue.seqCode

  residue = Residue(residue_type, number=residue_number)
  residue.name = residue_type

  if residue_number is not None:
      residue.name += str(residue_number)
  else:
      residue.name += str(i)
            
  ## create atoms

  for x in ccp_residue.atoms:
      residue.add_atom(make_isd_atom(x))

  ## link atoms

  cuts = connectivity.cuts

  for link in connectivity.links:

      name_1, name_2 = link.atom_1, link.atom_2

      if name_1 not in residue or name_2 not in residue:
          continue

      atom_1 = residue[name_1]
      atom_2 = residue[name_2]

      if link in cuts:
          atom_2.bind(atom_1)
      else:
          atom_2.link(atom_1)

  return residue

def create_isd_polymer(ccp_chain):
  """Descrn: Make an ISD polymer given a CCPN Chain object
     Inputs: ccp.molecule.MolSystem.Chain
     Output: ISD polymer
  """

  from Polymer_isd import Polymer
  from Connectivity import load_connectivity

  polymer = Polymer()

  connectivity = load_connectivity()
  
  isd_residues = [make_isd_residue(x, connectivity.getResidue(x.ccpCode)) for x in ccp_chain.sortedResidues()]

  ## link residues

  previous = isd_residues[0]

  for current in isd_residues[1:]:

      ## TODO: hard coded 'N' and 'C' for in and out going atoms

      current['N'].link(previous['C'])
      
      previous = current
      
  polymer[:] = isd_residues[:]
      
  polymer.define_groups()
  polymer.set_root(polymer.get_root())
  polymer.set_dofs()

  return polymer

def get_ccpn_constraint_list(ccpn_project, keys):
  """Descrn: Fetch a CCPN constraint list using list of object keys
     Inputs: List of keys (Words or Ints)
     Output: ccp.nmr.NmrConstraint.AbstractConstraintList
  """

  constraintList  = None
  storeSerial, serial = keys

  constraintStore = ccpn_project.findFirstNmrConstraintStore(serial=storeSerial)
  
  if not constraintStore:
      raise ValueError('No NMR constraint store with serial "%d" in CCPN project' % storeSerial)
  else:
      constraintList = constraintStore.findFirstConstraintList(serial=serial)

  if constraintList is None:
      raise ValueError('No constraint list with serial "%d" in store "%d"' % (serial,storeSerial))

  return constraintList

def get_ccpn_constraint_lists(ccpn_project, restraintsNames):
  """Descrn: Fetch a CCPN constraint list from a CCPN project using list 
             of restraint type and corresponding object keys
     Inputs: Implementation.Project, ARIA restraint names (what object is this?)
     Output: List if ccp.nmr.NmrConstraint.AbstractConstraintLists
  """

  constraintLists = {}
  for restraintType, constraintListKeys in restraintsNames.items():
    constraintLists[restraintType] = []

    for constraintListKey in constraintListKeys:
      keys = getKeysFromString(constraintListKey)
      constraintList = get_ccpn_constraint_list(ccpn_project, keys)

      if constraintList:
        constraintLists[restraintType].append(constraintList)
      
  return constraintLists

def make_struct_dict_from_polymer(polymer):

  d = {}

  seg_id = polymer.name

  altLoc = ' '
  insertCode = ' '

  for residue in polymer:

    d_residue = {}

    res_type = residue.type

    for name, atom in residue.items():

      d_atom = {'resName': res_type.capitalize(),
                'x': atom.x[0],
                'y': atom.x[1],
                'z': atom.x[2],
                'segId': seg_id,
                'insertCode': insertCode}

      d_residue[(name, altLoc)] = d_atom

    d[residue.number] = d_residue

  return {seg_id: d} 

def export_structure(mol_system, polymer, model_number):

  from ccpnmr.analysis.core.StructureIo import makeStructureEnsemble

  struct_dict = make_struct_dict_from_polymer(polymer)

  struct_dict = {model_number: struct_dict}

  try:
    structure = makeStructureEnsemble(struct_dict, mol_system) #, doWarnings=False)
  except:
    print struct_dict
    raise
  
  return structure

def remove_duplicates(l):

    d = {}

    for x in l:
        if not x in d:
            d[x] = 1

    return d.keys()

class CCPNReader:

    def __init__(self, filename=None, first_residue_number=1,
                 decompose_restraints=0, debug=False, project=None):

        from data import DataSet
        from Connectivity import load_connectivity

        self.debug = debug

        if project is not None:
          self.project = project
        else:
          self.project = self.open_ccpn_project(filename)

        self.project_filename = filename
        self.first_residue_number = first_residue_number

        self.connectivity = load_connectivity()

        self.decompose = decompose_restraints
        self.data_set = DataSet()
    
    def open_ccpn_project(self, filename):

        from memops.general.Io import loadProject
        import os

        filename = os.path.expanduser(filename)

        try:
            ccpn_project = loadProject(filename)

        except Exception, msg:
            raise IOError, 'Could not access CCPN project "%s". Error message was: %s' % (filename, msg)

        return ccpn_project

    def reset(self):

        import data

        sequence = self.data_set.sequence

        self.data_set = data.DataSet()
        self.data_set.sequence = sequence

    def resolve_dihedral_name(self, atoms, res_type):

        names = [a['name'] for a in atoms]

##         try:
##             res_type = self.sequence[atoms[1]['resid']]
                
##         except IndexError:
##             print 'Sequence not set or residue number overflow in atoms', atoms
##             return ''

        for dihedral in self.connectivity[res_type].dihedrals.values():

            keys = [k for k in dihedral.keys() if 'atom' in k]
            keys.sort()

            atom_names = []

            for k in keys:
                name = dihedral[k]
                if name[-1] in ('-', '+'):
                    name = name[:-1]

                atom_names.append(name)

            if atom_names == names:
                return dihedral['name']

        return None
                        
    def check_constraint_item(self, item, n_resonances=2):

        if len(item.resonances) <> n_resonances:
            raise ValueError, 'should be exactly %d resonances' % n_resonances

    def get_atom_list(self, resonance):

        if resonance.resonanceSet is None:
            return []

        l = [list(a.atoms) for a in resonance.resonanceSet.atomSets]

        return reduce(lambda a, b: a+b, l)

    def remove_duplicate_restraints(self, l):

        d = {}

        keys = []

        for x in l:

            key = tuple(x.contributions)

            if not key in d:
                d[key] = x
                keys.append(key)

            else:

                if self.debug:
                    print 'Discarding duplicate', x

        return [d[x] for x in keys]

    def get_volume(self, constraint, number):
        """
        attempts to find the volume of a ccpn constraint.
        """

        ## check whether we have reference peak

        value = None

        if constraint.peaks:

            peaks = constraint.peaks

            ## if we have more than one peak, take the first one.

            if len(peaks) > 1:
                if self.debug:
                    print 'Constraint %s: Multiple reference peaks, using first one.' % number
                    
            peak = list(peaks)[0]

            volume = peak.findFirstPeakIntensity(intensityType='volume')
            height = peak.findFirstPeakIntensity(intensityType='height')

            value = volume or height

            value = value.value
            
        if not value or not constraint.peaks:
            value = constraint.origData

            if value is not None and self.debug:
                print 'NOE or distance constraint %s: Using volume stored in original data.' % number

        if value is None:
            value = constraint.targetValue ** (-6.)
            if self.debug:
                print 'NOE or distance constraint %s: Using r^(-6) as volume.' % number

        return value

    def get_pair_contributions(self, constraint):

        contribs = []

        for item in constraint.items:

            self.check_constraint_item(item)

            resonances = list(item.resonances)

            atoms1 = self.get_atom_list(resonances[0])
            atoms2 = self.get_atom_list(resonances[1])

            if not atoms1 or not atoms2:
                continue

            for i in range(len(atoms1)):

                name1 = atoms1[i].name
                number1 = atoms1[i].residue.seqCode - \
                          self.first_residue_number

                for j in range(len(atoms2)):

                    name2 = atoms2[j].name
                    number2 = atoms2[j].residue.seqCode - \
                              self.first_residue_number

                    contribs.append(((number1, name1), (number2, name2)))

        return contribs
    
    def get_quad_contributions(self, constraint):

        if len(constraint.resonances) <> 4:

            print '4 resonances expected for DihedralConstraint, %d found' % len(constraint.resonances)
            return None, []
              
        atom_sets = [self.get_atom_list(R) for R in constraint.resonances]

        has_atoms = reduce(lambda a,b: a or b, atom_sets)

        if not has_atoms:
            return None, []

        contribs = []

        for atom_set in atom_sets:

            if not atom_set:
              return None, []

            name = atom_set[0].name
            number = atom_set[0].residue.seqCode - \
                   self.first_residue_number

            contribs.append((number, name))
            
        return [a[0] for a in atom_sets], contribs
    
    def convert_hbond_restraint_list(self, constraint_list):

        import data

        restraints = []

        for constraint in constraint_list.constraints:

            if self.debug:
              print 'HBond restraint %d:'% constraint.serial

            contribs = self.get_pair_contributions(constraint)

            if not contribs:
                if self.debug:
                    print 'No contributions found.'
                    
                continue

            contribs.sort()

            r_isd = data.Restraint()
            
            if constraint.targetValue is not None:
              r_isd.distance = constraint.targetValue

            else:
              
              if constraint.lowerLimit is not None and constraint.upperLimit is not None:

                r_isd.distance = 0.5*abs(constraint.upperLimit+constraint.lowerLimit)

                if not self.debug:
                  print 'HBond restraint %d:'% constraint.serial

                print 'No distance found. Using (upper_bound+lower_bound)/2 as estimate.'

              else:

                ## could display warning message

                pass

            r_isd.contributions = tuple(contribs)
##             r_isd.volume = self.get_volume(constraint, restraint_number)
            r_isd.upper = constraint.upperLimit
            r_isd.lower = constraint.lowerLimit
            
            r_isd.number = constraint.serial

            restraints.append(r_isd)
            
            if self.debug:
              print_distance_restraint(r_isd)

        r_list = data.DistanceList()
        r_list.restraints = restraints

        return r_list
 
    def convert_dihedral_restraint_list(self, constraint_list):
        
        import data
        from ccp.lib.MoleculeQuery import getAtomsTorsion

        restraints = []

        for constraint in constraint_list.constraints:

            if self.debug:
              print 'Dihedral restraint %d:' % constraint.serial

            atoms, contribs = self.get_quad_contributions(constraint)

            if not contribs:
                print 'No contributions found.'
                continue

            _atoms = [{'name': x.name, 'resid': x.residue.seqCode} for x in atoms]

            dihedral_name = self.resolve_dihedral_name(_atoms, atoms[1].residue.ccpCode)

            if not dihedral_name:
              print 'No torsion angle definition found for atoms %s' % str(tuple([(a.residue.seqCode, a.name) for a in atoms]))
              continue

##             torsion_obj = getAtomsTorsion(atoms)

##             if not torsion_obj:
##               print 'No torsion angle definition found for atoms %s' % str(tuple([(a.residue.seqCode, a.name) for a in atoms]))
##               continue

            res1 = atoms[1].residue.seqCode-self.first_residue_number
            res2 = atoms[2].residue.seqCode-self.first_residue_number

            if res1 <> res2:
              print 'Unable to define residue number unambiguously. Residue number for atoms 1 and 2: %d, %s.' % (res1, res2)
              continue

            res_number = res1
            
            if len(constraint.items) > 1:
              if self.debug:
                print 'Dihedral angle restraint %d has more than one DihedralConstraintItem. Using first one.' % constraint.serial

            try:
              item = constraint.items[0]
            except:
              print 'Skipping dihedral angle restraint #%d: No contributions found.' % constraint_index
              continue

            r_isd = data.TorsionAngleMeasurement()
            
            if item.targetValue is not None:
              r_isd.value = item.targetValue
              r_isd.error = item.error

            else:
              if item.lowerLimit is not None and item.upperLimit is not None:
                
                r_isd.value = 0.5*abs(item.upperLimit+item.lowerLimit)
                r_isd.error = 0.5*abs(item.upperLimit-item.lowerLimit)

                if not self.debug:
                  print 'Dihedral restraint %d:' % constraint.serial
                              
                print 'No target value found. Using (upper_bound+lower_bound)/2 as estimate.'

              else:

                ## could display warning here
                
                pass
              
##             r_isd.name = torsion_obj.name.lower()
            r_isd.name = str(dihedral_name)
            r_isd.residue_number = res_number
            r_isd.number = constraint.serial
            r_isd.contributions = (tuple(contribs),)

            restraints.append(r_isd)

            if self.debug:
              print 'Torsion angle %d: name=%s, residue number=%s, atoms=%s, value=%s' % \
                    (constraint.serial, r_isd.name, str(r_isd.residue_number), str(contribs), str(r_isd.value))

        r_list = data.TorsionAngleList()
        r_list.restraints = restraints

        return r_list

    def convert_jcoupling_restraint_list(self, constraint_list):
        
        import data
        from ccp.lib.MoleculeQuery import getAtomsTorsion
        from ccpnmr.analysis.ConstraintBasic import getConstraintAtoms

        restraints = []

        ## get list of jcouplings

        for constraint in constraint_list.constraints:

            atoms = getConstraintAtoms(constraint)

            contribs = self.get_quad_contributions(constraint)

            if not contribs:
                if self.debug:
                    print 'No contributions for restraint', constraint.serial
                    
                continue

            contribs = remove_duplicates(contribs)
            contribs.sort()

            r_isd = data.TorsionAngleMeasurement()
            r_isd.value = 0.
            r_isd.error = 0.
            r_isd.name = 'asdas'
            r_isd.redidue_number = 0
            
            r_isd.number = constraint.serial
            r_isd.contributions = tuple(contribs)
            r_isd.coupling = self.get_volume(constraint, restraint_number)

            restraints.append(r_isd)

        r_list = data.JCouplingList()
        r_list.restraints = restraints

        return r_list

    def convert_rdc_restraint_list(self, constraint_list):

        import data

        restraints = []

        for constraint in constraint_list.constraints:

            if constraint.targetValue is None and constraint.lowerLimit is None and \
               constraint.upperLimit is None:

              print 'No target value or bounds found for RDC', constraint.serial
              continue

            contribs = self.get_pair_contributions(constraint)

            if not contribs:
              print 'No contributions for RDC', constraint.serial
                    
              continue

            contribs = remove_duplicates(contribs)
            contribs.sort()

            r_isd = data.RDCMeasurement()
            r_isd.number = constraint.serial
            r_isd.contributions = tuple(contribs)

            if constraint.targetValue is not None:
              r_isd.coupling = constraint.targetValue

            else:

              print 'RDC restraint %d: No target value found, using (upper_bound+lower_bound)/2.' % constraint.serial

              r_isd.coupling = (constraint.upperLimit+constraint.lowerLimit)/2.

            restraints.append(r_isd)

        r_list = data.RDCList()
        r_list.restraints = restraints

        return r_list

    def convert_distance_restraint_list(self, constraint_list):

        import data

        restraints = []

        for constraint in constraint_list.constraints:
            
            if self.debug:
              print 'Distance restraint %d:'% constraint.serial
              
            contribs = self.get_pair_contributions(constraint)

            if not contribs:
                if self.debug:
                    print 'No contributions found.'
                    
                continue

            contribs = remove_duplicates(contribs)
            contribs.sort()

            r_isd = data.Restraint()
            r_isd.contributions = tuple(contribs)

            if constraint.targetValue is not None:
              r_isd.distance = constraint.targetValue

            else:
              
              if constraint.lowerLimit is not None and constraint.upperLimit is not None:

                r_isd.distance = 0.5*abs(constraint.upperLimit+constraint.lowerLimit)

                if not self.debug:
                  print 'HBond restraint %d:'% ccpn_restraint_number

                print 'No distance found. Using (upper_bound+lower_bound)/2 as estimate.'

              else:

                ## could display warning message

                pass

            r_isd.volume = self.get_volume(constraint, constraint.serial)
            r_isd.upper = constraint.upperLimit
            r_isd.lower = constraint.lowerLimit
            r_isd.number = constraint.serial

            restraints.append(r_isd)

        r_list = data.RestraintList()
        r_list.restraints = self.remove_duplicate_restraints(restraints)

        return r_list

    def _convert_restraint_list(self, constraint_list):

        key = constraint_list.className

        if key == 'DistanceConstraintList':
          
            restraint_list = self.convert_distance_restraint_list(constraint_list)
            
            if self.decompose:

              import data, TBLReader
              
              D = TBLReader.decompose_restraints(restraint_list.restraints)

              d = {}
                
              for key, restaints in d.items():
                R = data.RestraintList()
                R.restraints = restraints
                d[key] = R
                
            else:
                d = {None: restraint_list}

        elif key == 'HBondConstraintList':
            d = {None: self.convert_hbond_restraint_list(constraint_list)}

        elif key == 'RdcConstraintList':
            d = {None: self.convert_rdc_restraint_list(constraint_list)}

        elif key == 'DihedralConstraintList':
            d = {None: self.convert_dihedral_restraint_list(constraint_list)}

        elif key == 'JCouplingConstraintList':
            d = {None: self.convert_jcoupling_restraint_list(constraint_list)}

        else:
            print 'Unknown constrain list type "%s".' % key
            return None

        return d

    def write_xml(self, filename):

        from isdxml import ISDXMLPickler
        import os

        filename = os.path.expanduser(filename)

        pickler = ISDXMLPickler()
        pickler.dump(self.data_set, filename)

    def add_sequence(self, sequence, name):

        from data import Sequence

        seq = Sequence(name=name)

        for i in range(len(sequence)):
            seq.addResidue(i+self.first_residue_number, sequence[i])

        self.data_set.sequence = seq

    def read_sequence(self, key):

        keys = getKeysFromString(key)

        ccp_residues = get_ccpn_chain(self.project, keys).sortedResidues()

        sequence = [r.ccpCode for r in ccp_residues]

        self.first_residue_number = ccp_residues[0].seqCode
        
        self.add_sequence(sequence, key)

        d = {}

        for r in ccp_residues:
          d[r.seqCode] = r.ccpCode

        self.sequence = d

        return sequence

    def read_polymer(self, key):
        """
        creates and returns ISD polymer object
        """

        keys = getKeysFromString(key)

        ccp_chain = get_ccpn_chain(self.project, keys)

        polymer = create_isd_polymer(ccp_chain)

        return polymer

    def read_constraint_lists(self, keys):

        import data

        comment_template = 'CCPN constraint list. ISD identifier=%s'

        isd_restraint_lists = []

        for name in keys:

            if self.debug:
              print 'Reading restraint list %s' % name

            key = getKeysFromString(name)

            CL = get_ccpn_constraint_list(self.project, key)

            isd_restraints = self._convert_restraint_list(CL)

            if isd_restraints is None:
                continue

            l = []

            for label, r_list in isd_restraints.items():

                if not r_list.restraints:

                  if self.debug:
                    print 'Restraint list with key %s is empty.' % key
                    continue

                if label is not None:
                    r_list.key = '%s_%s' % (name, label)

                else:
                    r_list.key = name

                r_list.comment = comment_template % name

                self.data_set.add(r_list)

                isd_restraint_lists.append(r_list)

        return isd_restraint_lists

    def print_project_info(self):

        from cPickle import loads

        restraint_types = {'DistanceConstraintList': 'Distance/NOE',
                           'HBondConstraintList': 'H-Bond',
                           'RdcConstraintList': 'RDC',
                           'DihedralConstraintList': 'Dihedral angle',
                           'JCouplingConstraintList': 'J coupling'}

        ## sequences

        print '\nContents of CCPN project %s:\n' % self.project_filename

        print 'ISD projects:'

        template = ' - Key: %s (created on %s)'

        found = False

        for nmr_project in self.project.nmrProjects:

          app_data = nmr_project.findAllApplicationData(application=get_isd_name(),
                                                        keyword=key_project_settings)

          if not app_data:
            continue

          print '\nNmrProject: %s\n' % nmr_project.name

          for x in app_data:
            d = loads(x.value)

            date = d.get('date','N/A')
            
            print template % (repr(d['key']), date)

          found = True
          
        if found:
          print
        else:
          print '\n   None.\n'

        print 'Chains:\n'

        for mol_system in self.project.molSystems:
            for chain in mol_system.chains:

                key = getObjectKeyString(chain)

                print ' - ISD key: %s (CCPN details: %s)' % (repr(key), chain.details)

        ## constraint lists

        constraint_lists = []

        print '\nNMR projects:\n'

        for nmr_project in self.project.nmrProjects:
          key = getObjectKeyString(nmr_project)
          print ' - ISD key: %s (CCPN name: %s)' % (repr(key), nmr_project.name)
          constraint_lists += [x.sortedConstraintLists() for x in nmr_project.nmrConstraintStores]

        if constraint_lists:
          constraint_lists = reduce(lambda a,b: a+b, constraint_lists)

        print '\nConstraint lists: type (#restraints), ISD key (CCPN name):\n'

        for x in constraint_lists:

            key = getObjectKeyString(x)
            _type = restraint_types.get(x.className, 'Unknown')

            print ' - Type: %s (%.4d), ISD key: %s (CCPN name: %s)' % \
                  (_type, len(x.constraints), repr(key), x.name)

        if not constraint_lists:
          print '   None.'

        print

    def new_app_data_object(self, dest, src, obj_key, key=None):

      from cPickle import dumps
      from memops.api.Implementation import AppDataString
      from isd import VERSION_STRING
      import time

      if key is not None:

        d = {'key': key,
             'object': src,
             'date': time.ctime(),
             'isd_version': VERSION_STRING}

        string_val = dumps(d)
        
      else:
        string_val = dumps(src)
      
      app_data = AppDataString(value=string_val,
                               application=get_isd_name(),
                               keyword=obj_key)

      dest.addApplicationData(app_data)

      return app_data

    def get_app_data_object(self, app_data):

      import cPickle

      return cPickle.loads(app_data.value)['object']

    def find_app_data(self, app_data, key):

      import cPickle

      for x in app_data:

        if type(x.value) <> str:
          continue

        try:
          d = cPickle.loads(x.value)
        except cPickle.UnpicklingError:
          continue
          
        if d['key'] == key:
          return x

      return None

    def export_project_settings(self, settings, nmr_project_name, key):

      print 'Exporting ISD project to CCPN NmrProject %s using key "%s" ...' % \
            (nmr_project_name, key)

      nmr_project = self.get_nmr_project(nmr_project_name, True)

      ## make sure application data with key "key" does not exist

      app_data = self.find_app_data(nmr_project.applicationData, key)

      if app_data is None:
        app_data = self.new_app_data_object(nmr_project, settings, key_project_settings, key)

      else:

##         for x in nmr_project.applicationData:
##           nmr_project.removeApplicationData(x)

        print 'Exising project settings with key "%s" will be overwritten.' % key

        nmr_project.removeApplicationData(app_data)
        app_data = self.new_app_data_object(nmr_project, settings, key_project_settings, key)

    def get_nmr_project(self, nmr_project_name, create_new=False):

      nmr_project = list(self.project.findAllNmrProjects(name=nmr_project_name))

      if len(nmr_project) > 1:
          raise KeyError, 'Key %s for CCPN NmrProject is not unique.' % \
                nmr_project_name

      elif not nmr_project:

        if create_new:
          print 'CCPN NmrProject "%s" not found. Creating new project.' % nmr_project_name

          from ccp.api.nmr.Nmr import NmrProject
          
          nmr_project = NmrProject(self.project, name=nmr_project_name)
          
        else:
          raise KeyError, 'CCPN NmrProject with key %s not found.' % \
                nmr_project_name

      else:
        nmr_project = nmr_project[0]

      return nmr_project

    def read_project_settings(self, nmr_project_name, key):

      from cPickle import loads

      print 'Reading ISD project settings from CCPN NmrProject %s using key %s  ...' % \
            (nmr_project_name, key)

      nmr_project = self.get_nmr_project(nmr_project_name)

      app_data = nmr_project.findAllApplicationData(application=get_isd_name(),
                                                    keyword=key_project_settings)

      if not app_data:
        raise KeyError, 'CCPN NmrProject %s contains no ISD simulation settings with key %s.' % \
              (nmr_project_name, key)
      
      ad = self.find_app_data(app_data, key)

      if ad is None:
        raise KeyError, 'CCPN NmrProject %s contains no ISD simulation settings with key %s.' % \
              (nmr_project_name, key)

      return self.get_app_data_object(ad)

    def export_ensemble(self, mol_system_name, ensemble, index_list, polymer):

        from ccpnmr.analysis.core.StructureIo import makeStructureEnsemble

        ## get molSystem

        mol_systems = list(self.project.findAllMolSystems(code=mol_system_name))

        if not mol_systems:

            print 'CCPN MolSystem %s not found. Creating new one.' % \
                  mol_system_name

            from ccp.api.molecule.MolSystem import MolSystem

            mol_system = MolSystem(self.project, code=mol_system_name)

        elif len(mol_systems) <> 1:
            raise KeyError, 'More than one CCPN MolSystem with name %s found.' % \
                  mol_system_name

        else:
            mol_system = mol_systems[0]

        structures = {}

        for i in range(len(index_list)):

            index = index_list[i]

            print index

            if self.debug:
                print 'Exporting structure index=%d' % index

            state = ensemble[index]

            polymer.set_torsions(state.torsion_angles, 1)

            ## create CCPN MolStructure

            structures[i] = make_struct_dict_from_polymer(polymer)
            
        ensemble = makeStructureEnsemble(structures, mol_system)

    def save(self):
      if self.debug:
        print 'Saving CCPN project ...'
        
      self.project.saveModified()

    def set_sequence(self, d):
      self.sequence = d

if __name__ == '__main__':

  import sys, os

  modules = os.environ['ISD_ROOT'] + '/src/py'

  if not modules in sys.path:
    sys.path.insert(0, modules)

  from setup import create_simulation_from_project as F
  from utils import Load

  ccpn_project = '/home/wr222/projects/aria2.1/examples/werner/ccpn/hrdc_ccpn_project.xml'
  ccpn_project = '/home/wr222/test_isd/tudor.xml'
  ccpn_project = '/home/wr222/projects/data/tudor/ccpn/tudor.xml'

  mol_system_key = 'hrdc|A'
  mol_system_key = 'Molecularsystem|A'
  constraint_keys = ('tudor|1|1',)
  constraint_keys = ('tudor|1|1', 'tudor|1|2', 'tudor|1|3', 'tudor|1|4')
  constraint_keys = ('tudor|1|1', 'tudor|1|2', 'tudor|1|2', 'tudor|1|4', 'tudor|1|5',
                     'tudor|1|6', 'tudor|1|7', 'tudor|1|8')

##   simulation = F('~/simulations/tudor/floating_test/floating_test')
##  simulation = Load('/tmp/sim')
##   E = Load('~/simulations/tudor/floating_test/floating_test_0')

  C = CCPNReader(ccpn_project, decompose_restraints=False, debug=True)
##   C.export_simulation(simulation, 'tudor')
##   sim = C.import_simulation('tudor')
##   C.export_ensemble('tudor', E, [1,], simulation.posterior.get_polymer())

##       C.print_project_info()

##       polymer = C.read_sequence(mol_system_key)

  CL = C.read_constraint_lists(constraint_keys)

##       C.write_xml('filename')

  from memops.api.Implementation import Project

  p = Project('test')
