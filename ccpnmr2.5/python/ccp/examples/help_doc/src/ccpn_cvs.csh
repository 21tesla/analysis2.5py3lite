# Or whatever your SourceForge username is.
setenv SF_USER "penkett"

# Or wherever you want to put the CCPN repository.
setenv TOP_DIR `pwd`

setenv CCPNMR_TOP_DIR "${TOP_DIR}/ccpn"
setenv PYTHONPATH ".:${CCPNMR_TOP_DIR}/python"

cd $TOP_DIR
cvs -d :pserver:${SF_USER}@ccpn.cvs.sourceforge.net:/cvsroot/ccpn \
    co -r HEAD ccpn

cd ccpn/python/memops/scripts_v2/
python makePython.py
