include ../../environment.txt
# above gives CC, LINK, MAKE, CC_FLAGS, LINK_FLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB

GLOBAL_DIR = ../../memops/global
GLOBAL_INCLUDE_FLAGS = -I$(GLOBAL_DIR)

GLOBAL_PY_UTIL_OBJS = \
	$(GLOBAL_DIR)/utility.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/python_util.$(OBJ_SUFFIX)

PY_BACUS_OBJS = \
	bacus.$(OBJ_SUFFIX) \
	py_bacus.$(OBJ_SUFFIX)

all: Bacus.$(DYLIB_SUFFIX)

Bacus.$(DYLIB_SUFFIX): global_py_util_objects $(PY_BACUS_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_BACUS_OBJS) $(MATH_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

py_bacus.$(OBJ_SUFFIX): py_bacus.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

.c.$(OBJ_SUFFIX):; $(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(CO_NAME)

global_py_util_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_util_objects

clean:
	$(RM) *.$(OBJ_SUFFIX)

realclean: clean
	$(RM) *.$(DYLIB_SUFFIX)
