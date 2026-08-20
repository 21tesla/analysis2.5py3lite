include environment.txt

all: compile

compile: global structure analysis clouds dynamics bayes

global:
	cd memops/global && $(MAKE)

structure:
	cd ccp/structure && $(MAKE)

analysis:
	cd ccpnmr/analysis && $(MAKE)

clouds:
	cd ccpnmr/clouds && $(MAKE)

dynamics:
	cd ccpnmr/dynamics && $(MAKE)

bayes:
	cd other/cambridge/bayes && $(MAKE)

clean: global_clean structure_clean analysis_clean clouds_clean dynamics_clean bayes_clean

global_clean:
	cd memops/global && $(MAKE) clean

structure_clean:
	cd ccp/structure && $(MAKE) clean

analysis_clean:
	cd ccpnmr/analysis && $(MAKE) clean

clouds_clean:
	cd ccpnmr/clouds && $(MAKE) clean

dynamics_clean:
	cd ccpnmr/dynamics && $(MAKE) clean

bayes_clean:
	cd other/cambridge/bayes && $(MAKE) clean

realclean: global_realclean structure_realclean analysis_realclean clouds_realclean dynamics_realclean bayes_realclean

global_realclean:
	cd memops/global && $(MAKE) realclean

structure_realclean:
	cd ccp/structure && $(MAKE) realclean

analysis_realclean:
	cd ccpnmr/analysis && $(MAKE) realclean

clouds_realclean:
	cd ccpnmr/clouds && $(MAKE) realclean

dynamics_realclean:
	cd ccpnmr/dynamics && $(MAKE) realclean

bayes_realclean:
	cd other/cambridge/bayes && $(MAKE) realclean

links: memops_links ccp_links ccpnmr_links other_links

memops_links:
	cd ../python/memops/c && $(LINK_LIBRARIES)

ccp_links:
	cd ../python/ccp/c && $(LINK_LIBRARIES)

ccpnmr_links:
	cd ../python/ccpnmr/c && $(LINK_LIBRARIES)

other_links:
	cd ../python/cambridge/c && $(LINK_LIBRARIES)

copies: memops_copies ccp_copies ccpnmr_copies other_copies

memops_copies:
	cd ../python/memops/c && $(COPY_LIBRARIES)

ccp_copies:
	cd ../python/ccp/c && $(COPY_LIBRARIES)

ccpnmr_copies:
	cd ../python/ccpnmr/c && $(COPY_LIBRARIES)

other_copies:
	cd ../python/cambridge/c && $(COPY_LIBRARIES)

cleanpyd: memops_pyd ccp_pyd ccpnmr_pyd other_pyd

memops_pyd:
	cd ../python/memops/c && $(CLEAN_PYD)

ccp_pyd:
	cd ../python/ccp/c && $(CLEAN_PYD)

ccpnmr_pyd:
	cd ../python/ccpnmr/c && $(CLEAN_PYD)

other_pyd:
	cd ../python/cambridge/c && $(CLEAN_PYD)

cleanso: memops_so ccp_so ccpnmr_so other_so

memops_so:
	cd ../python/memops/c && $(CLEAN_SO)

ccp_so:
	cd ../python/ccp/c && $(CLEAN_SO)

ccpnmr_so:
	cd ../python/ccpnmr/c && $(CLEAN_SO)

other_so:
	cd ../python/cambridge/c && $(CLEAN_SO)

#test:
#	$(PYTHON_DIR)/bin/python 

