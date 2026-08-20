"""

Wrapper around Gary's Marvin code to fit into formatConverter structure

"""
# TODO: needs better code for the handling of heavyAtomShifts both from and to
import os, string, re

# Import general functions
from ccp.format.marvin.generalIO import MarvinGenericFile
from memops.universal.Io import splitPath
from memops.universal.Util import returnInt

#####################
# Class definitions #
#####################

class MarvinPeakFile(MarvinGenericFile):

  """
  Information on file level
  """

  def initialize(self,assignTagSep = ' '):

    self.peaks = []
    self.peakIdPatt = re.compile('3d[c|n]')

    self.peakById = {}

  def setSpectrumInfo(self,specName,ndim):

    self.specNames = [specName]
    self.numDims = [ndim]




  def read(self,verbose = 0):
    if verbose == 1:
      print "Reading %s peak list %s" % (self.format,self.name)



    # Initialize
    # ----------

    (path,specName) = splitPath(self.name)
    self.setSpectrumInfo(specName,3)

    fin = open(self.name)
    reader = Marvin_reader()
    reader.read(fin)

    handler = Print_assignments(returnValue = True)

    # TODO: This is not correct - should be able to import the whole peak list if required!

    for peak in reader.records:
      if 'peakAssign' in peak:
        for assignment in peak['peakAssign']:
          likelihood=assignment['likelihood'][0].as_float()
          # hack indeting problems and unneeded
          if 1:
          #if likelihood > probability:
            #0       1                 2                3               4                        5
            #peak_id,from_proton_shift,from_heavy_shift,to_proton_shift,from_assignment.shift_id,from_assignment.residue_id,
            #6                            7                         8                      9
            #from_assignment.residue_type,from_assignment.atom_name,to_assignment.shift_id,to_assignment.residue_id,
            #10                         11                      12          13                            14
            #to_assignment.residue_type,to_assignment.atom_name,probability,from_assignment.residue_type, to_assignment.residue_type

            values = handler.handle_assignment(peak,assignment)
            peakId = returnInt(values[0])

            # Dangerous, works for most things but really needs a table driven approach

            heavy_name = None
            #convertAtomNameByType(self,values[13],values[7])


            light_name =  values[7]
#            if light_name== 'HN':
#                heavy_name='N'
#            elif light_name[-1].isdigit() or light_name[-1] == '#':
#                heavy_name = light_name[:-1]
#                heavy_name=heavy_name.replace('H','C')
#            else:
#                heavy_name=light_name.replace('H','C')

            if values[2] != None:
                assign = [self.getAssignString(values[5],'LIGHT:'+values[13]+':'+values[7]),
                          self.getAssignString(values[5],'HEAVY:'+values[13]+':'+values[7]),
                          self.getAssignString(values[9],'LIGHT:'+values[14]+':'+values[11])]
                ppm = values[1:4]
                label = [values[6],values[6],values[10]]
            else:
                assign = [self.getAssignString(values[5],'LIGHT:'+values[13]+':'+values[7]),
                          self.getAssignString(values[9],'LIGHT:'+values[14]+':'+values[11])]
                ppm = [values[1],values[3]]
                label = [values[6],values[10]]

            #print values[0],assign, values[5],values[7],values[9],values[11],heavy_name


            if peakId not in self.peakById:
              newPeak = MarvinPeak(peakId,assign,ppm,label,likelihood)
              print (peakId,assign,ppm,label,likelihood)
              self.peaks.append(newPeak)
              self.peakById[peakId] = self.peaks[-1]
            else:
              self.peakById[peakId].addAssignment(assign,label,likelihood)

            #print assign
    fin.close()

  def getAssignString(self,seqCode,atom):

    return "%s.%s.%s" % (self.defaultMolCode,seqCode,atom)

class MarvinPeak:

  def __init__(self,num,assign,ppm,label,probability):

    self.num = num
    self.ppm = ppm
    self.assignList = [assign]
    self.label = [label]
    self.probabilities = [probability]

  def addAssignment(self,assign,label,probability):

    self.assignList.append(assign)
    self.label.append(label)
    self.probabilities.append(probability)

"""

Code by Gary Thompson

"""

import sys
import re
from UserDict import DictMixin
# todo: rename values

class Record(object,DictMixin):
    def __init__(self,type,values,leaf=False):
        self.leaf=leaf
        self.type=type
        self.values=values
        self.children = []

    def dump(self,level=1):
        result= []
        result.append('%s : %s' % (self.type,   self.values))
        for elem in self.children:
            indent = '\t'*level
            for line in  elem.dump(level=level+1).split('\n'):
                result.append('%s%s' % (indent,line))
        return '\n'.join(result)

    def as_float(self):
        return float(self.values)

    def as_string(self):
        return `self.values`

    def is_container(self):
        return len(self.children)

    def __str__(self):
        return self.dump()

#    def __setitem__(self,key,value):
#        raise NotImplementedError('__setitem__ is not implimented this is a read only dictionary')
#
#    def __delitem__(self,key):
#        raise NotImplementedError('__delitem__ is not implimented this is a read only dictionary')

    def keys(self):
        key_set=set()
        for elem in self.children:
            key_set.add(elem.type)

        return list(key_set)


    def __getitem__(self,key):
        result = []
        for elem in self.children:
            if elem.type == key:
                result.append(elem)

        if len(result) != 0:
            return result
        else:
            raise KeyError(key)

class Marvin_reader(object):
    def __init__(self,handler=None,prefix='3d[c|n]'):
        self.handler=handler
        self.records=[]
        self.prefix_re=re.compile(prefix)

    def read(self,file):


        for record in self.pass_data(file):

            self.records.append(record)



    # todo: need dump handler as default
    def dump(self, handler=None,sort=True,ostream=sys.stdout):

        sort_list = []
        if sort:
            for record in self.records:
                peak_id=self.prefix_re.sub('',record.values)
                sort_list.append((int(peak_id),record))
            sort_list.sort()
        self.records = [pair[1] for pair in sort_list]

        if handler == None:
            print >> ostream, '! dumping marvin peak records'
            print >> ostream, '!---------------------------'
            print >> ostream, ''
            print >> ostream, '!%d marvin peak records' %  len(self.records)
            print >> ostream, ''
            for record in self.records:
                print >> ostream, record
                print >> ostream, ''
        else:

            for record in self.records:
                handler.handle_peak(record)




    def pass_data(self,file):
      for line in file:
           bare_line = self.clear_comments(line).strip()

           if len(bare_line) > 0:
               fields = bare_line.split()
               yield self.process_record(type=fields[0],fields=' '.join(fields[1:]),file=file)



    def process_record(self,type,fields, file):
        record = Record(type,fields)
        for line in file:
            bare_line = self.clear_comments(line).strip()

            if len(bare_line) == 0:
                continue
            fields = bare_line.split()

            if bare_line[0] == '-':
                record.children.append(Record(type=fields[0][1:],values=' '.join(fields[1:]),leaf=True))

            elif fields[0] == 'end':
                break
            else:
                child= self.process_record(type=fields[0],fields=' '.join(fields[1:]),file=file)
                record.children.append(child)
        return record

    def clear_comments(self, line):
        pat=re.compile('([^!])*!.*$')

        match  = re.match(pat, line)
        if match:
            line = match.group(1)
            if line:
                line = line.strip()
            else:
                line = ''

        return line

class Peak_handler(object):
    def handle_peak(self,record):
        pass

class Assignment_handler(object):
    def handle_assignment(self,peak,record):
        pass

class Peak_assignment_adapter(Peak_handler):
    def __init__(self,assignment_handler):
        self.assignment_handler=assignment_handler

    def handle_peak(self,peak):
        if 'peakAssign' in peak:
            for assignment in peak['peakAssign']:
                self.assignment_handler.handle_assignment(peak,assignment)


class Assignment_probability_filter(Assignment_handler):
    def __init__(self,handler,probability=0.9):
        probability=float(probability)
        assert probability >= 0.0
        assert probability <= 1.0
        self.probability=probability
        self.handler=handler

    def handle_assignment(self,peak,assignment):
         likelihood=assignment['likelihood'][0].as_float()
         if likelihood > self.probability:
             self.handler.handle_assignment(peak,assignment)

def Dump_assignment(Assignment_handler):
    def __init__(self,ostream=sys.__stdout__):
        self.ostream=ostream

    def handle_assignment(record):
        print >> self.ostream, peak_id,from_proton_shift,from_heavy_shift,to_proton_shift,

        assignment= assignment.values
        (from_assignment,to_assignment) = self.values_to_assignments(assignment)

        print >> self.ostream, from_assignment.shift_id,from_assignment.residue_id,from_assignment.residue_type,from_assignment.atom_name,
        print >> self.ostream, to_assignment.shift_id,to_assignment.residue_id,to_assignment.residue_type,to_assignment.atom_name,
        print >> self.ostream, merit




class Dump_peak_handler(Peak_handler):
    def handle_peak(self,record):
        print record.dump()

class Dump_filter(Peak_handler):
    def __init__(self,sub_handler):
        self.sub_handler=sub_handler

    def check_peak_ok(self,record):
        pass

    def handle_peak(self,record):
        if self.check_peak_ok(record):
            self.sub_handler.handle_peak(record)

class Assigned_peak_filter(Dump_filter):
    def __init__(self,sub_handler):
        super(Assigned_peak_filter,self).__init__(sub_handler)

    def check_peak_ok(self,record):
        return 'peakAssign' in record




class Print_assignments(Assignment_handler):
    class Axis_assignment(object):
        def __init__(self,shift_id,residue_id,residue_type,atom_name):
            self.shift_id=shift_id
            self.residue_id=residue_id
            self.residue_type=residue_type
            self.atom_name=atom_name

    class Assignment_pair(object):
        def __init__(self,from_assignment,to_assignment):
           self.from_assignment=from_assignment
           self.to_assignment=to_assignment

    def __init__(self,prefix='3d[c|n]',ostream=sys.stdout,returnValue = False):
        self.prefix_re=re.compile(prefix)
        self.ostream=ostream
        self.returnValue = returnValue

    def decompose_assignment_string(self,string):

        string=string.replace('_from(',' ')
        string=string.replace('_to(',' ')
        string=string.replace(')','')
        string=string.replace('_',' ')


        (shift_id,residue_id,residue_type,atom_name) = string.split(' ')
        shift_id=self.prefix_re.sub('',shift_id)

        return self.Axis_assignment (shift_id,residue_id,residue_type,atom_name)

    def values_to_assignments(self,string):
        (dummy_peak_assignment,from_proton,to_proton) = string.split()
        dummy_peak_assignment=dummy_peak_assignment #hack to avoid usused error

        return (self.decompose_assignment_string(from_proton),self.decompose_assignment_string(to_proton))

    def handle_assignment(self,peak,assignment):
        peak_id = peak.values
        peak_id=self.prefix_re.sub('',peak_id)

        from_proton_shift = peak['fromProtonShift'][0].as_float()
        try:
            from_heavy_shift = peak['fromHeavyatomShift'][0].as_float()
        except:
            from_heavy_shift =  None

        to_proton_shift = peak['toProtonShift'][0].as_float()

        assignment_name= assignment.values
        (from_assignment,to_assignment) = self.values_to_assignments(assignment_name)

        probability = assignment['likelihood'][0].as_float()


        if not self.returnValue:

          print >> self.ostream, peak_id,from_proton_shift,from_heavy_shift,to_proton_shift,
          print >> self.ostream, from_assignment.shift_id,from_assignment.residue_id,from_assignment.residue_type,from_assignment.atom_name,
          print >> self.ostream, to_assignment.shift_id,to_assignment.residue_id,to_assignment.residue_type,to_assignment.atom_name,
          print >> self.ostream, probability

        else:
          #print from_assignment.residue_id,to_assignment.residue_id
          return (peak_id,from_proton_shift,from_heavy_shift,to_proton_shift,
                  from_assignment.shift_id,from_assignment.residue_id,from_assignment.residue_type,from_assignment.atom_name,
                  to_assignment.shift_id,to_assignment.residue_id,to_assignment.residue_type,to_assignment.atom_name,
                  probability,from_assignment.residue_type,to_assignment.residue_type)
