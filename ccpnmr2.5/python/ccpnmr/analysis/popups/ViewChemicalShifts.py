
"""
======================COPYRIGHT/LICENSE START==========================

ViewChemicalShifts.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""

from ccpnmr.analysis.core.AssignmentBasic import getResonanceName, getShiftLists
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames, getResidueCode

from memops.gui.ButtonList import UtilityButtonList, ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Label import Label
from memops.gui.Frame import Frame
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledCanvas import ScrolledCanvas

colorDict = {'H':'#a0d0a0','N':'#a0a0d0',
             'C':'#d0d0a0','F':'#d0d0d0',
             'P':'#d0a0d0','S':'#d0c090',}


def ChemicalShiftsTableMacro(argServer, shiftList=None):

  popup = ViewChemicalShiftsPopup(argServer.parent)
  
  popup.open()
    
    
class ViewChemicalShiftsPopup(BasePopup):
  """
  **A Table of Chemical Shifts for Export**
  
  This section is designed to make a layout of a table for chemical shifts  on a
  per-residue basis which may them be exported as either PostScript, for
  printing and graphical manipulation, or as plain text for import into other
  software or computer scripts. 

  The user chooses the molecular chain (which sequence) and the shift list to
  use at the top of the popup, together with a few other options that control
  how things are rendered. Then buttons are toggled to select which kinds of
  atom will be displayed in aligned columns; other kinds will simply be listed
  to the right of the columns. Thus for example if the shift list does not
  contain any carbonyl resonances in a protein chain then the user may toggle
  the empty "C" column off.

  Once the desired layout is achieved the user then uses the [Export PostScript]
  or [Export Text] buttons to write the data into a file of the appropriate
  type. The user will be presented wit ha file browser to specify the location
  and the name of the file to be saved. It should be noted that although the
  graphical display in the popup itself is somewhat limited, e.g. the gaps and
  spacing doesn't always look perfect, the PostScript version that is exported
  is significantly neater.

  **Caveats & Tips**

  If you need a chemical shift list represented in a particular format, specific
  for a particular external NMR program then you should use the FormatConverter
  software.

  Chemical shifts may also be exported from any table in Analysis that contains
  such data by clicking the right mouse button over the table and selecting the
  export option. 

  """

  def __init__(self, parent, *args, **kw):

    self.shiftList  = None
    self.font       = 'Helvetica 10'
    self.boldFont   = 'Helvetica 10 bold'
    self.symbolFont = 'Symbol 8'
    self.smallFont  = 'Helvetica 8'
    self.chain      = None
    self.textOut    = ''
    
    BasePopup.__init__(self, parent=parent, title='Chart : Chemical Shifts Table')

  def body(self, guiFrame):
  
      
    row = 0
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(None,6)
    
    label = Label(frame, text='Chain:', grid=(0,0))
    tipText = 'Selects which molecular chain to show residues and chemical shift values for'
    self.chainPulldown = PulldownList(frame, callback=self.changeChain,
                                      grid=(0,1), tipText=tipText)
    
    label = Label(frame, text='  Shift List:', grid=(0,2))
    tipText = 'Selects which shift list is used to derive the displayed chemical shift values'
    self.shiftListPulldown = PulldownList(frame, callback=self.changeShiftList,
                                          grid=(0,3), tipText=tipText)

    label = Label(frame, text=' List all shifts:', grid=(0,4))
    tipText = 'Sets whether to display all the chemical shifts for residues or just for the nominated atom types in columns'
    self.otherShiftsSelect = CheckButton(frame, callback=self.draw,
                                         grid=(0,5), tipText=tipText)

    utilButtons = UtilityButtonList(frame, helpUrl=self.help_url,
                                    grid=(0,7))

    row += 1
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(None,6)
    
    label = Label(frame, text=' 1-letter codes:', grid=(0,0))
    tipText = 'Whether to use 1-letter residue codes in the table, or otherwise Ccp/three-letter codes'
    self.oneLetterSelect = CheckButton(frame, callback=self.draw,
                                       grid=(0,1), selected=False, tipText=tipText)
    
    precisions = [0.1,0.01,0.001]
    texts = [str(t) for t in precisions]
    label = Label(frame, text='  1H precision:', grid=(0,2))
    tipText = 'Specifies how many decimal places to use when displaying 1H chemical shift values'
    self.protonPrecisionSelect = PulldownList(frame, texts=texts, objects=precisions,
                                              callback=self.draw, index=1, grid=(0,3),
                                              tipText=tipText)

    label = Label(frame, text='  Other precision:')
    label.grid(row=0, column=4, sticky='w')
    tipText = 'Specifies how many decimal places to use when displaying chemical shift values for isotopes other than 1H'
    self.otherPrecisionSelect = PulldownList(frame, texts=texts, objects=precisions,
                                             callback=self.draw, index=1, grid=(0,5),
                                             tipText=tipText)
 
    row  += 1
    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(None,1)
    
    label = Label(frame, text='Column\nAtoms:', grid=(0,0))
    tipText = 'Selects which kinds of atoms are displayed in aligned columns, or otherwise displayed at the end of the residue row (if "List all shifts" is set)'
    self.optSelector=PartitionedSelector(frame, self.toggleOpt, tipText=tipText,
                                         maxRowObjects=10, grid=(0,1), sticky='ew')
    options = ['H','N','C','CA','CB','CG']
    self.optSelector.update(objects=options, labels=options,
                            selected=['H','N','CA'])
    
    row  += 1
    guiFrame.expandGrid(row, 0)
    self.canvasFrame = ScrolledCanvas(guiFrame,relief = 'groove', width=650,
                                      borderwidth=2, resizeCallback=None,
                                      grid=(row,0), padx=1, pady=1 )
    self.canvas = self.canvasFrame.canvas
    #self.canvas.bind('<Button-1>', self.toggleResidue)
  
    row += 1
    tipTexts = ['Output information from the table as PostScript file, for printing etc.',
                'Output information from the table as a whitespace separated plain text file']
    commands = [self.makePostScript,self.exportText]
    texts    = ['Export PostScript','Export Text'] 
    buttonList = ButtonList(guiFrame, commands=commands, texts=texts,
                            grid=(row,0), tipTexts=tipTexts)
    
    chains = self.getChains()
    if len(chains) > 1:
      self.chain = chains[1]
    else:  
      self.chain = None
    
    self.updateShiftLists()
    self.updateChains()
    self.otherShiftsSelect.set(True)
    self.update()
    
    for func in ('__init__','delete'):
      self.registerNotify(self.updateChains,'ccp.molecule.MolSystem.Chain',func)
    for func in ('__init__','delete'):
      self.registerNotify(self.updateShiftLists,'ccp.nmr.Nmr.ShiftList',func)

  
  def changeShiftList(self, shiftList):
      
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      self.update()
  
  def updateShiftLists(self, *opt):
  
    names = []
    index = 0
  
    shiftLists = getShiftLists(self.nmrProject)
    
    shiftList  = self.shiftList
    
    if shiftLists:
      names =  ['%s [%d]' % (sl.name or '<No name>', sl.serial) for sl in shiftLists]
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
        
      index = shiftLists.index(shiftList)  
    
    if self.shiftList is not shiftList:
      self.shiftList = shiftList
           
    self.shiftListPulldown.setup(names,shiftLists,index)
    
 
  def getChains(self):

    chains = [None,]
    for molSystem in self.project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        if len(chain.residues) > 1:
          chains.append(chain)
    
    return chains    

  def changeChain(self, chain):
        
    if chain is not self.chain:
      self.chain = chain
      self.update()

  def updateChains(self, *opt):
  
    names = []
    index = -1
  
    chains = self.getChains()
    chain  = self.chain
    
    if len(chains) > 1:
      names = ['None',] + ['%s:%s' % (ch.molSystem.code,ch.code) for ch in chains[1:]]
      if chain not in chains:
        chain = chains[1]
      
      index = chains.index(chain)  
       
    else:
      chain = None
       
    self.chainPulldown.setup(names, chains, index)
  
  def destroy(self):

    for func in ('__init__','delete'):
      self.unregisterNotify(self.updateChains,'ccp.molecule.MolSystem.Chain',func)
    for func in ('__init__','delete'):
      self.unregisterNotify(self.updateShiftLists,'ccp.nmr.Nmr.ShiftList',func)
  
    BasePopup.destroy(self)

  def exportText(self, *event):
  
    from memops.gui.FileSelect import FileType
    from memops.gui.FileSelectPopup import FileSelectPopup

    if self.textOut:
      fileTypes = [  FileType('Text', ['*.txt']), FileType('All', ['*'])]
      fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
                 title = 'Save table as plain text', dismiss_text = 'Cancel',
                 selected_file_must_exist = False)

      fileName = fileSelectPopup.getFile()
    
      if fileName:
        file = open(fileName, 'w')
        file.write(self.textOut)
  
  def toggleOpt(self, selected):
  
    self.draw()
  
  def makePostScript(self):
  
    self.canvasFrame.printCanvas()
    
  def update(self):
  
    colors  = []
    selected = set()
    atomNames = set()
 
    if self.chain:
      for residue in self.chain.sortedResidues():
        for atom in residue.atoms:
          chemAtom = atom.chemAtom
          
          if colorDict.get(chemAtom.elementSymbol) is None:
            continue
        
          if chemAtom.waterExchangeable:
            continue
        
          atomNames.add(atom.name[:2])
 
        molType = residue.molResidue.molType
        if molType == 'protein':
          selected.add('H')
          selected.add('N')
          selected.add('CA')
          selected.add('C')
          
        elif molType == 'DNA':
          selected.add("C1'")
          
        elif molType == 'RNA':
          selected.add("C1'")
          
        elif molType == 'carbohydrate':
          selected.add("C1")
 
    else:
      for spinSystem in self.shiftList.nmrProject.resonanceGroups:
        if not spinSystem.residue:
          for resonance in spinSystem.resonances:
            for name in resonance.assignNames:
               atomNames.add(name)
 
    options = list(atomNames)
    
    molType = 'protein'
    if self.chain:
      molType = self.chain.molecule.molType
    options = greekSortAtomNames(options, molType=molType)
    
    if 'H' in options:
      options.remove('H')
      options = ['H',] + options
    colors  = [colorDict.get(n[0]) for n in options]
 
    if options and not selected:
      selected = [options[0],]
 
    self.optSelector.update(objects=options,
                            labels=options,
                            colors=colors,
                            selected=list(selected))
 
    self.draw()
 
  def draw(self, *opt):
  
    if not self.shiftList:
      return
  
    nmrProject = self.shiftList.nmrProject
    shiftList  = self.shiftList
    font       = self.font
    bfont      = self.boldFont
    symbolFont = self.symbolFont
    sFont      = self.smallFont
    bbox       = self.canvas.bbox
    doOthers   = self.otherShiftsSelect.get()
  
    spc = 4
    gap = 14
    x   = gap
    y   = gap
    ct  = self.canvas.create_text
    cl  = self.canvas.create_line
    cc  = self.canvas.coords
    
    self.canvas.delete('all')
  
    ssDict = {}
  
    formatDict = {0.1:'%.1f',0.01:'%.2f',0.001:'%.3f',}
  
    protonFormat = formatDict[self.protonPrecisionSelect.getObject()]
    otherFormat  = formatDict[self.otherPrecisionSelect.getObject()]
    
    uSpinSystems = []
    chains       = set()
    molSystems   = set()
    
    for spinSystem in nmrProject.resonanceGroups:
 
      residue = spinSystem.residue
      if residue:
        ssDict[residue] = ssDict.get(residue, []) + [spinSystem,]
      else:
        uSpinSystems.append( (spinSystem.serial, spinSystem) )
  
    uSpinSystems.sort()

    
    commonAtoms = self.optSelector.getSelected()
    N = len(commonAtoms)
    
    chain = self.chain
      
    if chain:
      spinSystems = []
      for residue in chain.sortedResidues():
        spinSystems0 = ssDict.get(residue, [])
 
        for spinSystem in spinSystems0:
          if spinSystem and spinSystem.resonances:
            spinSystems.append([residue, spinSystem])
    else:
      spinSystems = uSpinSystems
    
   
    strings = []
    
    doOneLetter = self.oneLetterSelect.get()
    
    if spinSystems:
     
      x = gap
      y += gap
 
      numItems    = []
      codeItems   = []
      commonItems = []
      otherItems  = []
 
      numWidth  = 0
      codeWidth = 0
      commonWidths = [0] * N
      commonCounts = [0] * N
      
      for residue, spinSystem in spinSystems:
        
        if type(residue) is type(1):
          seqNum  = '{%d}' % residue
          
          if doOneLetter:
            ccpCode = '-'
          else:
            ccpCode = spinSystem.ccpCode or ''
          
        else:
          if doOneLetter:
            ccpCode = residue.chemCompVar.chemComp.code1Letter
          
          else:
            ccpCode = getResidueCode(residue)
 
          seqNum = str(residue.seqCode)
        
        subStrings = []
        subStrings.append(seqNum)
        subStrings.append(ccpCode)
        
        item = ct(x,y,text=seqNum,font=font,anchor='se')
        box  = bbox(item)
        iWidth = box[2]-box[0]
        numWidth = max(numWidth,iWidth)
        numItems.append(item)

        item = ct(x,y,text=ccpCode,font=font,anchor='sw')
        box  = bbox(item)
        iWidth = box[2]-box[0]
        codeWidth = max(codeWidth,iWidth)
        codeItems.append(item)
 
        commonShifts, commonElements, otherShifts = self.getShiftData(spinSystem, shiftList, commonAtoms)
 
        items = []
        for i in range(N):
          values = commonShifts[i]
          element =  commonElements[i]
          
          if element == 'H':
            shiftFormat = protonFormat
          else:
            shiftFormat = otherFormat
          
          subItems = []
          for value in values:
            text = shiftFormat % value
                     
            if text:
              item = ct(x,y,text=text,font=font,anchor='se')
              box  = bbox(item)
              iWidth = box[2]-box[0]
              commonWidths[i] = max(commonWidths[i],iWidth)
              commonCounts[i] += 1
              subItems.append(item)           
 
          subStrings.append(','.join([shiftFormat % v for v in values]) or '-')
  
          items.append(subItems)
         
        commonItems.append(items)
        
        if doOthers:
          items0 = []
          i = 0
          I = len(otherShifts)
          for atomLabel, element, value in otherShifts:
          
            label = atomLabel
            if label[0] == '?':
              label = label[3:-1]
            
            if element == 'H':
              shiftFormat = protonFormat
            else:
              shiftFormat = otherFormat
            
            subStrings.append('%6s:%-4s' % (shiftFormat % value,label))
            i += 1
 
            atoms = atomLabel.split('|')
 
            items = []
            j = 0
            for atom in atoms:
              text = element
              if j > 0:
                text = '/' + text
 
              item = ct(x,y,text=text,font=font,anchor='sw')
              box  = bbox(item)
              iWidth = box[2]-box[0]-3
              items.append((iWidth, item, 0))
 
              p = len(element)
              if len(atom) > p:
                letter = atom[p]
                if letter not in ('0123456789'):
                  item = ct(x,y,text=letter.lower(),font=symbolFont,anchor='sw')
                  box  = bbox(item)
                  iWidth = box[2]-box[0] -2
                  items.append((iWidth, item, -4))
                  p += 1
 
              text = atom[p:]
              if text:
                item = ct(x,y,text=text,font=sFont,anchor='sw')
                box  = bbox(item)
                iWidth = box[2]-box[0]-2
                items.append((iWidth, item, -4))
              j += 1
 
            text = ' ' + shiftFormat % value
            if i != I:
              text += ','
 
            item = ct(x,y,text=text,font=font,anchor='sw')
            box  = bbox(item)
            iWidth = box[2]-box[0] - 4
            items.append((iWidth, item,0))
 
            items0.append(items)
        
          otherItems.append(items0)
          
          
        strings.append(subStrings)
 
      y0 = y
      x = x0 = gap + numWidth + codeWidth + spc + spc
      for i in range(N):
        
        if not commonCounts[i]:
          continue
        
        x += commonWidths[i] + spc + spc
        
        
        element = commonElements[i]
        
        iWidth = 0
        text = commonAtoms[i][len(element):].lower()
        if text:
          item = ct(x,y,text=text,font=symbolFont,anchor='se')
          box  = bbox(item)
          iWidth = box[2]-box[0]-2

        ct(x-iWidth,y,text=element,font=font,anchor='se')
      
      y += gap
        
      for i in range(len(numItems)):
        x = gap + numWidth + spc
        
        cc(numItems[i],x,y)
      
        x += spc
        
        cc(codeItems[i],x,y)

        x += codeWidth
        
        x1 = x + spc
        
        yM = y
        for j in range(N):
          if not commonCounts[j]:
            continue
            
          x += commonWidths[j] + spc + spc
        
          items = commonItems[i][j]
          
          yB = y-gap
          for item in items:
            yB += gap
            cc(item, x, yB)
            yM = max(yB,yM)
        
        x += gap

        if doOthers:
          x += spc
          x3 = x
 
          for items in otherItems[i]:
            if x > 550:
              x = x3
              y += gap
 
            for iWidth, item, dy in items:
              cc(item, x, y+dy)
              x += iWidth
      
        y = max(y,yM)
        y += gap
      
      x = x0
      for i in range(N):
        if not commonCounts[i]:
          continue
          
        x += commonWidths[i] + spc + spc
          
        cl(x+8,y0,x+8,y-gap,width=0.3,fill='#808080')
        
      cl(x1,y0,x1,y-gap,width=0.3,fill='#808080')
      cl(0,0,550,0,width=0.3,fill='#FFFFFF')
        
        
      y += gap
    
    textWidths = {} 
    for subStrings in strings:
      for i, text in enumerate(subStrings):
        if text and len(text) > textWidths.get(i, 0):
          textWidths[i] = len(text)
        else:
          textWidths[i] = 0
        
    formats = {}
    for i in textWidths.keys():
      formats[i] = ' %%%ds' % max(6,textWidths[i])

    textOut = '!' 
    textOut += ' '*(max(6,textWidths.get(0,6))+max(6,textWidths.get(1,6))+1)
    
    i = 2
    for atom in commonAtoms:
      if i in formats:
        textOut += formats[i] % atom
        i += 1
    textOut += '\n'
    
    for subStrings in strings:
      i = 0
      for text in subStrings:
        textOut += formats[i] % text
        i += 1
      
      textOut += '\n'
       
    self.textOut = textOut
         
  def getShiftData(self, spinSystem, shiftList, commonAtoms= ('H','N','CA','CB')):

    commonShifts     = []
    commonResonances = {}
    commonElements   = []
 
    for atomName in commonAtoms:
      elements = set()
      resonances = []
 
      for resonance in spinSystem.resonances:
        resonanceSet = resonance.resonanceSet

        if resonanceSet:
          for atomSet in resonanceSet.atomSets:
            for atom in atomSet.atoms:
              if atomName == atom.name[:2]:
                 resonances.append(resonance)
                 break
                 
            else:
              continue
            break       
 
        else:
          for assignName in resonance.assignNames:
            if atomName == assignName[:2]:
              resonances.append(resonance)
              break

      shiftValues = []
      for resonance in resonances:
        isotope = resonance.isotope
        if isotope:
          elements.add(isotope.chemElement.symbol)
        commonResonances[resonance] = True
        shift = resonance.findFirstShift(parentList=shiftList)
        if shift:
          shiftValues.append(shift.value)

      if not elements:
        element = atomName[0]
      else:
        element = elements.pop()

      commonElements.append(element)
      commonShifts.append(shiftValues)

    otherShifts = []
    for resonance in spinSystem.resonances:
      if not commonResonances.get(resonance):
        shift = resonance.findFirstShift(parentList=shiftList)
 
        if shift:
          isotope = resonance.isotope
          
          if isotope:
            element = isotope.chemElement.symbol
          else:
            element = '??'
 
          if resonance.assignNames or resonance.resonanceSet:
            name = getResonanceName(resonance)
          else:
            name = '??[%d]' % resonance.serial

          otherShifts.append((name, element, shift.value))

    molType = 'protein'
    if spinSystem.residue:
      molType = spinSystem.residue.molResidue.molType

    otherShifts = greekSortAtomNames(otherShifts, molType)

    return commonShifts, commonElements, otherShifts
