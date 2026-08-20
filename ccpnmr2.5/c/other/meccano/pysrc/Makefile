include ../../../environment.txt
# above gives CC, CC_FLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB

GLOBAL_DIR = ../../../memops/global
GLOBAL_INCLUDE_FLAGS = -I$(GLOBAL_DIR)

MECCANO_MAIN_DIR = ../meccano2_stat_ramaDB_fwd
MECCANO_INC_DIR = ../inc
MECCANO_SRC_DIR = ../src
MECCANO_INCLUDE_FLAGS = -I$(MECCANO_INC_DIR)

GLOBAL_PY_UTIL_OBJS = \
	$(GLOBAL_DIR)/python_util.o \
	$(GLOBAL_DIR)/utility.o
 
MECCANO_OBJS = \
	$(MECCANO_SRC_DIR)/gradient.o \
	$(MECCANO_SRC_DIR)/myDAT.o \
	$(MECCANO_SRC_DIR)/myGEOMETRY.o \
	$(MECCANO_SRC_DIR)/myHBSC.o \
	$(MECCANO_SRC_DIR)/myJUNC.o \
	$(MECCANO_SRC_DIR)/myMINIMISATION.o \
	$(MECCANO_SRC_DIR)/myPPGEO.o \
	$(MECCANO_SRC_DIR)/myRAMACHANDRAN.o \
	$(MECCANO_SRC_DIR)/myRDC.o \
	$(MECCANO_SRC_DIR)/mySTRUCT.o \
	$(MECCANO_MAIN_DIR)/meccano2_stat_ramaDB_fwd.o
 
PY_MECCANO_OBJS = \
	py_meccano.o

all: Meccano.so

Meccano.so: $(PY_MECCANO_OBJS) global_py_util_objects meccano_objects
	$(CC) $(SHARED_FLAGS) -o $@ $(GSL_LIB_FLAGS) $(PYTHON_LIB_FLAGS) \
		$(PY_MECCANO_OBJS) $(GLOBAL_PY_UTIL_OBJS) $(MECCANO_OBJS) \
		$(GSL_LIB) $(PYTHON_LIB) -lm

py_meccano.o: py_meccano.c
	$(CC) -c $(GLOBAL_INCLUDE_FLAGS) $(GSL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(MECCANO_INCLUDE_FLAGS) $(CC_FLAGS) $<

.c.o:;  $(CC) -c $(GLOBAL_INCLUDE_FLAGS) $(MECCANO_INCLUDE_FLAGS) $(CC_FLAGS) $<
 
global_py_util_objects:
	cd $(GLOBAL_DIR); $(MAKE) py_util_objects

meccano_objects:
	cd $(MECCANO_MAIN_DIR); $(MAKE) meccano_objects

clean:
	rm *.o

realclean:
	rm *.o
	rm *.so

