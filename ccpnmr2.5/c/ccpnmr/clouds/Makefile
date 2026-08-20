include ../../environment.txt
# above gives CC, LINK, MAKE, CC_FLAGS, LINK_FLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB

GLOBAL_DIR = ../../memops/global
GLOBAL_INCLUDE_FLAGS = -I$(GLOBAL_DIR)

GLOBAL_PY_UTIL_OBJS = \
	$(GLOBAL_DIR)/utility.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/python_util.$(OBJ_SUFFIX)

GLOBAL_RANDOM_OBJS = \
	$(GLOBAL_DIR)/random.$(OBJ_SUFFIX)

GLOBAL_DIAG_OBJS = \
	$(GLOBAL_DIR)/diag_dbl.$(OBJ_SUFFIX)

DYNAMICS_OBJS = \
	atom_coord.$(OBJ_SUFFIX) \
	atom_coord_list.$(OBJ_SUFFIX) \
	dist_constraint.$(OBJ_SUFFIX) \
	dist_constraint_list.$(OBJ_SUFFIX) \
	dist_force.$(OBJ_SUFFIX) \
	dynamics.$(OBJ_SUFFIX)

PY_DYNAMICS_OBJS = \
	$(DYNAMICS_OBJS) \
	py_atom_coord.$(OBJ_SUFFIX) \
	py_atom_coord_list.$(OBJ_SUFFIX) \
	py_dist_constraint.$(OBJ_SUFFIX) \
	py_dist_constraint_list.$(OBJ_SUFFIX) \
	py_dist_force.$(OBJ_SUFFIX) \
	py_dynamics.$(OBJ_SUFFIX)

PY_ATOM_COORD_OBJS = \
	atom_coord.$(OBJ_SUFFIX) \
	py_atom_coord.$(OBJ_SUFFIX)

PY_ATOM_COORD_LIST_OBJS = \
	$(PY_ATOM_COORD_OBJS) \
	atom_coord_list.$(OBJ_SUFFIX) \
	py_atom_coord_list.$(OBJ_SUFFIX)

PY_DIST_CONSTRAINT_OBJS = \
	dist_constraint.$(OBJ_SUFFIX) \
	py_dist_constraint.$(OBJ_SUFFIX)

PY_DIST_CONSTRAINT_LIST_OBJS = \
	$(PY_DIST_CONSTRAINT_OBJS) \
	dist_constraint_list.$(OBJ_SUFFIX) \
	py_dist_constraint_list.$(OBJ_SUFFIX)

PY_DIST_FORCE_OBJS = \
	dist_force.$(OBJ_SUFFIX) \
	py_dist_force.$(OBJ_SUFFIX)

PY_MIDGE_OBJS = \
	midge.$(OBJ_SUFFIX) \
	py_midge.$(OBJ_SUFFIX)

PY_BACUS_OBJS = \
	bacus.$(OBJ_SUFFIX) \
	py_bacus.$(OBJ_SUFFIX)

PY_CLOUD_UTIL_OBJS = \
	py_cloud_util.$(OBJ_SUFFIX)

all: AtomCoord.$(DYLIB_SUFFIX) AtomCoordList.$(DYLIB_SUFFIX) DistForce.$(DYLIB_SUFFIX) DistConstraint.$(DYLIB_SUFFIX) \
	DistConstraintList.$(DYLIB_SUFFIX) Dynamics.$(DYLIB_SUFFIX) Midge.$(DYLIB_SUFFIX) CloudUtil.$(DYLIB_SUFFIX) Bacus.$(DYLIB_SUFFIX)

AtomCoord.$(DYLIB_SUFFIX): global_py_util_objects $(PY_ATOM_COORD_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_ATOM_COORD_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

AtomCoordList.$(DYLIB_SUFFIX): global_py_util_objects $(PY_ATOM_COORD_LIST_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_ATOM_COORD_LIST_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DistConstraint.$(DYLIB_SUFFIX): global_py_util_objects $(PY_DIST_CONSTRAINT_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_DIST_CONSTRAINT_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DistConstraintList.$(DYLIB_SUFFIX): global_py_util_objects $(PY_DIST_CONSTRAINT_LIST_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_DIST_CONSTRAINT_LIST_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DistForce.$(DYLIB_SUFFIX): global_py_util_objects $(PY_DIST_FORCE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_DIST_FORCE_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

Dynamics.$(DYLIB_SUFFIX): global_random_objects global_py_util_objects $(PY_DYNAMICS_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_RANDOM_OBJS) $(GLOBAL_PY_UTIL_OBJS) \
		$(PY_DYNAMICS_OBJS) $(MATH_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

Midge.$(DYLIB_SUFFIX): global_diag_objects global_py_util_objects $(PY_MIDGE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_DIAG_OBJS) $(GLOBAL_PY_UTIL_OBJS) \
		$(PY_MIDGE_OBJS) $(MATH_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

Bacus.$(DYLIB_SUFFIX): global_py_util_objects $(PY_BACUS_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_BACUS_OBJS) $(MATH_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

CloudUtil.$(DYLIB_SUFFIX): $(PY_CLOUD_UTIL_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(PY_CLOUD_UTIL_OBJS) $(MATH_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

py_atom_coord.$(OBJ_SUFFIX): py_atom_coord.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_atom_coord_list.$(OBJ_SUFFIX): py_atom_coord_list.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_dist_constraint.$(OBJ_SUFFIX): py_dist_constraint.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_dist_constraint_list.$(OBJ_SUFFIX): py_dist_constraint_list.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_dist_force.$(OBJ_SUFFIX): py_dist_force.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_dynamics.$(OBJ_SUFFIX): py_dynamics.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_midge.$(OBJ_SUFFIX): py_midge.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_bacus.$(OBJ_SUFFIX): py_bacus.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_cloud_util.$(OBJ_SUFFIX): py_cloud_util.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

.c.$(OBJ_SUFFIX):; $(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(CO_NAME)
 
global_py_util_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_util_objects

global_diag_objects:
	cd $(GLOBAL_DIR) && $(MAKE) diag_objects

global_random_objects:
	cd $(GLOBAL_DIR) && $(MAKE) random_objects

clean:
	$(RM) *.$(OBJ_SUFFIX)

realclean: clean
	$(RM) *.$(DYLIB_SUFFIX)

