
import os

#ccpnProjectNames = ('AR3436A','CGR26A','CtR69A','ET109A','HR5537A','NeR103A',
#                    'PGR122A','VpR247','AtT13','HR6470A','HR6430A','HR5460A','OR36','OR135',
#                    'StT322','HR2876B','YR313A','MiR12','HR8254A')
ccpnProjectNames = ('HR6470A','HR6430A','HR5460A','OR36','OR135',
                    'StT322','HR2876B','YR313A','HR8254A','HR2876C')

#skipEntries = (
# 139, 197, 248, 280, 306, 292,  # Superseded 
 #135, 136, 207, 277, 279,      # Duplicate LYS 1 HZ records, FIXED
 #281,                          # incorrect methyl names, FIXED
#)


casdNmrDir = os.environ.get('CASD_HOME')
if not casdNmrDir:
  raise Exception("Environment variable CASD_HOME not set")
allDataDir = os.path.join(casdNmrDir, 'data')
topTmpDir = os.path.join(casdNmrDir, 'tmpdata')
if not os.path.exists(topTmpDir):
  os.makedirs(topTmpDir)

#casdResultsDir = os.path.join(casdNmrDataDir, 'results')
#casdInputDir = os.path.join(casdNmrDataDir, 'input')
#casdInputPdbDir = os.path.join(casdInputDir, 'pdb')
#casdInputCcpnDir = os.path.join(casdInputDir, 'ccpn')
#casdDownloadDir = os.path.join(casdInputDir, 'downloads')

pdbEndings = ('.pdb', '.pdb_v3.1')

