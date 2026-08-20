include ../../environment.txt
# above gives CC, LINK, MAKE, CC_FLAGS, LINK_FLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB

GLOBAL_DIR = ../../memops/global
GLOBAL_INCLUDE_FLAGS = -I$(GLOBAL_DIR)

GLOBAL_PY_UTIL_OBJS = \
	$(GLOBAL_DIR)/utility.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/python_util.$(OBJ_SUFFIX)

GLOBAL_RANDOM_OBJS = \
	$(GLOBAL_DIR)/random.$(OBJ_SUFFIX)

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

all: DyAtomCoord.$(DYLIB_SUFFIX) DyAtomCoordList.$(DYLIB_SUFFIX) DyDistForce.$(DYLIB_SUFFIX) DyDistConstraint.$(DYLIB_SUFFIX) \
	DyDistConstraintList.$(DYLIB_SUFFIX) DyDynamics.$(DYLIB_SUFFIX)

DyAtomCoord.$(DYLIB_SUFFIX): global_py_util_objects $(PY_ATOM_COORD_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_ATOM_COORD_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DyAtomCoordList.$(DYLIB_SUFFIX): global_py_util_objects $(PY_ATOM_COORD_LIST_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_ATOM_COORD_LIST_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DyDistConstraint.$(DYLIB_SUFFIX): global_py_util_objects $(PY_DIST_CONSTRAINT_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_DIST_CONSTRAINT_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DyDistConstraintList.$(DYLIB_SUFFIX): global_py_util_objects $(PY_DIST_CONSTRAINT_LIST_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_DIST_CONSTRAINT_LIST_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DyDistForce.$(DYLIB_SUFFIX): global_py_util_objects $(PY_DIST_FORCE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_DIST_FORCE_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

DyDynamics.$(DYLIB_SUFFIX): global_random_objects global_py_util_objects $(PY_DYNAMICS_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_RANDOM_OBJS) $(GLOBAL_PY_UTIL_OBJS) \
		$(PY_DYNAMICS_OBJS) $(MATH_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

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

.c.$(OBJ_SUFFIX):; $(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(CO_NAME)
 
global_py_util_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_util_objects

global_random_objects:
	cd $(GLOBAL_DIR) && $(MAKE) random_objects

clean:
	$(RM) *.$(OBJ_SUFFIX)

realclean: clean
	$(RM) *.$(DYLIB_SUFFIX)

