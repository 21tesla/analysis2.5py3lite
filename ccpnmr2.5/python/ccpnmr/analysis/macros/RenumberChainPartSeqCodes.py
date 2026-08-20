def renumberChainPartSeqCodes(argServer):
  
  chain = argServer.getChain()
  
  if not chain: # None avail or aborted
    msg = 'Stopping: Script cancelled or no chain available.'
    argServer.showWarning(msg)
    return

  currFirstSeqCode = argServer.askInteger('Enter starting seqCode to renumber:', 1)
  
  if currFirstSeqCode is None:
    msg = 'Stopping: Script cancelled.'
    argServer.showWarning(msg)
    return 
  
  newFirstSeqCode = argServer.askInteger('Enter what seqCode that should now be:', 1)
  
  if newFirstSeqCode is None:
    msg = 'Stopping: Script cancelled.'
    argServer.showWarning(msg)
    return 
  
  gap = newFirstSeqCode - currFirstSeqCode

  for residue in chain.residues:
    seqCode = residue.seqCode
    if seqCode >= currFirstSeqCode:
      residue.seqCode += gap

  argServer.showInfo('Chain renumbered successfully')
