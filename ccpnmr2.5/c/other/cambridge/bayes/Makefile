# ==================================================================
# = Makefile for the NMR BayeSys program =
# = Created by Daniel O'Donovan on 2008-09-22.
# = Copyright (c) 2008 University of Cambridge. All rights reserved.
# ==================================================================

include ../../../environment.txt
# above gives CC, LINK, MAKE, CC_FLAGS, LINKFLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB, OPENMP_*, MATH_LIB

GLOBAL_DIR = ../../../memops/global
GLOBAL_INCLUDE_FLAGS = -I$(GLOBAL_DIR)

GLOBAL_PY_UTIL_OBJS = \
	$(GLOBAL_DIR)/utility.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/python_util.$(OBJ_SUFFIX)

GLOBAL_PY_MEM_CACHE_OBJS = \
	$(GLOBAL_PY_UTIL_OBJS) \
	$(GLOBAL_DIR)/hash_list.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/hash_table.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/mem_cache.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/mutex.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_mem_cache.$(OBJ_SUFFIX)

GLOBAL_PY_BLOCK_FILE_OBJS = \
	$(GLOBAL_PY_MEM_CACHE_OBJS) \
	$(GLOBAL_DIR)/block_file.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/shape_file.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/int_array.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_block_file.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_shape_file.$(OBJ_SUFFIX)


LOCAL_FLAGS	= $(CC_FLAGS) $(OPENMP_FLAGS)
# LOCAL_FLAGS	= -g $(OPENMP_FLAGS)

BAYES_OBJ	=   random.$(OBJ_SUFFIX) hilbert.$(OBJ_SUFFIX) bayesys3.$(OBJ_SUFFIX)

BNMR_PY_OBJ	=	bayes_nmr.$(OBJ_SUFFIX) 
BNMR_OBJ	=	bayes_nmr.$(OBJ_SUFFIX) 
PP_OBJ		=	app.$(OBJ_SUFFIX) distribution.$(OBJ_SUFFIX) 

PY_OBJ		=	$(GLOBAL_PY_BLOCK_FILE_OBJS) py_bayes.$(OBJ_SUFFIX) $(BNMR_PY_OBJ)

LOCAL_INC		= $(PYTHON_INCLUDE_FLAGS) $(GLOBAL_INCLUDE_FLAGS)
LOCAL_LIB		= $(PYTHON_LIB) $(PLOT_LIB) $(OPENMP_LIB) $(MATH_LIB)

#PY_BAYES_OBJ	=   $(BAYES_OBJ) $(PY_OBJ) $(PP_OBJ) $(PLOT_OBJ)
PY_BAYES_OBJ	=   $(BNMR_OBJ) $(PP_OBJ) $(BAYES_OBJ) py_bayes.$(OBJ_SUFFIX)

SHARED_LIB		=	BayesPeakSeparator.$(DYLIB_SUFFIX)

#all: $(PY_OBJ) $(BAYES_OBJ) $(PP_OBJ) $(GLOBAL_PY_BLOCK_FILE_OBJS) $(SHARED_LIB)
all: $(SHARED_LIB)

# all files (general compile)
.c.$(OBJ_SUFFIX):;  $(CC) $(LOCAL_FLAGS) $(LOCAL_INC) $(CO_NAME)

$(SHARED_LIB): global_py_block_file_objects $(PY_BAYES_OBJ)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(GLOBAL_LIB_FLAGS) $(PYTHON_LIB_FLAGS) \
	$(GLOBAL_PY_BLOCK_FILE_OBJS) $(PY_BAYES_OBJ) $(PYTHON_LIB) $(LOCAL_LIB)

global_py_util_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_util_objects

global_py_mem_cache_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_mem_cache_objects

global_py_block_file_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_block_file_objects

clean:
	$(RM) *.$(OBJ_SUFFIX)

realclean: clean
	$(RM) *.$(DYLIB_SUFFIX)
