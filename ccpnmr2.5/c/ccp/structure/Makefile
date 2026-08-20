include ../../environment.txt
# above gives CC, LINK, MAKE, CC_FLAGS, LINKFLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB

GLOBAL_DIR = ../../memops/global
GLOBAL_INCLUDE_FLAGS = -I$(GLOBAL_DIR)

GLOBAL_COLOR_OBJS = \
	$(GLOBAL_DIR)/color.$(OBJ_SUFFIX)

GLOBAL_PY_UTIL_OBJS = \
	$(GLOBAL_DIR)/utility.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/python_util.$(OBJ_SUFFIX)
 
GLOBAL_SORTS_OBJS = \
	$(GLOBAL_DIR)/sorts.$(OBJ_SUFFIX)

GLOBAL_GEOMETRY_OBJS = \
	$(GLOBAL_DIR)/geometry.$(OBJ_SUFFIX)

GLOBAL_LINALG_OBJS = \
	$(GLOBAL_DIR)/eigenvalue.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/linalg.$(OBJ_SUFFIX)

GLOBAL_DRAWING_OBJS = \
	$(GLOBAL_COLOR_OBJS) \
	$(GLOBAL_DIR)/clipping.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/gl_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/pdf_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/ps_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/store_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/tk_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/w32func.$(OBJ_SUFFIX)

GLOBAL_PY_DRAWING_OBJS = \
	$(GLOBAL_DRAWING_OBJS) \
	$(GLOBAL_DIR)/py_draw_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_gl_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_pdf_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_ps_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_store_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_tk_handler.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_tk_util.$(OBJ_SUFFIX)

ALL_INCLUDE_FLAGS = $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) \
	$(X11_INCLUDE_FLAGS) \
	$(TCL_INCLUDE_FLAGS) $(TK_INCLUDE_FLAGS) $(GL_INCLUDE_FLAGS)

ALL_LIB_FLAGS = $(X11_LIB_FLAGS) $(TK_LIB_FLAGS) $(GL_LIB_FLAGS)

ALL_LIB = $(GL_LIB) $(TK_LIB) $(X11_LIB) $(MATH_LIB)

PY_STRUCT_ATOM_OBJS = \
	atom.$(OBJ_SUFFIX) \
	bond.$(OBJ_SUFFIX) \
	py_atom.$(OBJ_SUFFIX)

PY_STRUCT_BOND_OBJS = \
	$(PY_STRUCT_ATOM_OBJS) \
	py_bond.$(OBJ_SUFFIX)

PY_STRUCT_STRUCTURE_OBJS = \
	$(PY_STRUCT_BOND_OBJS) \
	structure.$(OBJ_SUFFIX) \
	py_structure.$(OBJ_SUFFIX)

PY_STRUCT_UTIL_OBJS = \
	struct_util.$(OBJ_SUFFIX) \
	py_struct_util.$(OBJ_SUFFIX)

all: StructAtom.$(DYLIB_SUFFIX) StructBond.$(DYLIB_SUFFIX) StructStructure.$(DYLIB_SUFFIX) StructUtil.$(DYLIB_SUFFIX)

StructAtom.$(DYLIB_SUFFIX): global_color_objects global_py_util_objects $(PY_STRUCT_ATOM_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(GLOBAL_LIB_FLAGS) $(PYTHON_LIB_FLAGS) \
			$(GLOBAL_COLOR_OBJS) $(GLOBAL_PY_UTIL_OBJS) $(PY_STRUCT_ATOM_OBJS) $(PYTHON_LIB)

StructBond.$(DYLIB_SUFFIX): global_color_objects global_py_util_objects $(PY_STRUCT_BOND_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(GLOBAL_LIB_FLAGS) $(PYTHON_LIB_FLAGS) \
			$(GLOBAL_COLOR_OBJS) $(GLOBAL_PY_UTIL_OBJS) $(PY_STRUCT_BOND_OBJS) $(PYTHON_LIB)

StructStructure.$(DYLIB_SUFFIX): global_geometry_objects global_sorts_objects \
			global_py_util_objects global_py_drawing_objects $(PY_STRUCT_STRUCTURE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(ALL_LIB_FLAGS) $(PYTHON_LIB_FLAGS) \
			$(GLOBAL_GEOMETRY_OBJS) $(GLOBAL_SORTS_OBJS) \
			$(GLOBAL_PY_UTIL_OBJS) $(GLOBAL_PY_DRAWING_OBJS) $(PY_STRUCT_STRUCTURE_OBJS) \
			$(PYTHON_LIB) $(ALL_LIB)

StructUtil.$(DYLIB_SUFFIX): global_geometry_objects global_linalg_objects global_py_util_objects $(PY_STRUCT_UTIL_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(GLOBAL_LIB_FLAGS) $(PYTHON_LIB_FLAGS) \
			$(GLOBAL_GEOMETRY_OBJS) $(GLOBAL_LINALG_OBJS) \
			$(GLOBAL_PY_UTIL_OBJS) $(PY_STRUCT_UTIL_OBJS) $(PYTHON_LIB)

py_atom.$(OBJ_SUFFIX): py_atom.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_bond.$(OBJ_SUFFIX): py_bond.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_structure.$(OBJ_SUFFIX): py_structure.c
	$(CC) $(CC_FLAGS) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_struct_util.$(OBJ_SUFFIX): py_struct_util.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

.c.$(OBJ_SUFFIX):;  $(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(CO_NAME)
 
global_py_util_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_util_objects

global_color_objects:
	cd $(GLOBAL_DIR) && $(MAKE) color_objects

global_sorts_objects:
	cd $(GLOBAL_DIR) && $(MAKE) sorts_objects

global_geometry_objects:
	cd $(GLOBAL_DIR) && $(MAKE) geometry_objects

global_linalg_objects:
	cd $(GLOBAL_DIR) && $(MAKE) linalg_objects

global_drawing_objects:
	cd $(GLOBAL_DIR) && $(MAKE) drawing_objects

global_py_drawing_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_drawing_objects

clean:
	$(RM) *.$(OBJ_SUFFIX)

realclean: clean
	$(RM) *.$(DYLIB_SUFFIX)

