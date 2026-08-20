def renumberChainSeqCodes(argServer):
  
  chain = argServer.getChain()
  
  if not chain: # None avail or aborted
    msg = 'Stopping: Script cancelled or no chain available.'
    argServer.showWarning(msg)
    return

  firstSeqCode = argServer.askInteger('Enter starting seqCode:', 1)
  
  if firstSeqCode is None:
    msg = 'Stopping: Script cancelled.'
    argServer.showWarning(msg )
    return 
  
  if firstSeqCode < 0 and firstSeqCode+len(chain.residues) >= 0:
    skipZeroSeqCode = argServer.askYesNo('Do you want to skip 0 in numbering?')
  else:
    skipZeroSeqCode = False # arbitrary
  
  from ccp.util.Molecule import renumberChainSeqCodes as renumber

  renumber(chain, firstSeqCode, skipZeroSeqCode)

  argServer.showInfo('Chain renumbered successfully')
