include ../../environment.txt
# above gives CC, LINK, MAKE, CC_FLAGS, LINK_FLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB

GLOBAL_DIR = ../../memops/global
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
 
GLOBAL_PY_STORE_FILE_OBJS = \
	$(GLOBAL_DIR)/store_file.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/py_store_file.$(OBJ_SUFFIX)
 
GLOBAL_CONTOUR_OBJS = \
	$(GLOBAL_DIR)/contourer.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/list.$(OBJ_SUFFIX)

GLOBAL_NONLINEAR_MODEL_OBJS = \
	$(GLOBAL_DIR)/nonlinear_model.$(OBJ_SUFFIX) \
	$(GLOBAL_DIR)/gauss_jordan.$(OBJ_SUFFIX)

GLOBAL_DRAWING_OBJS = \
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

PY_CONTOUR_FILE_OBJS = \
	contour_data.$(OBJ_SUFFIX) \
	contour_levels.$(OBJ_SUFFIX) \
	contour_style.$(OBJ_SUFFIX) \
	contour_file.$(OBJ_SUFFIX) \
	py_contour_levels.$(OBJ_SUFFIX) \
	py_contour_style.$(OBJ_SUFFIX) \
	py_contour_file.$(OBJ_SUFFIX)

PY_SLICE_FILE_OBJS = \
	slice_file.$(OBJ_SUFFIX) \
	py_slice_file.$(OBJ_SUFFIX)

PY_CONTOUR_STYLE_OBJS = \
	contour_style.$(OBJ_SUFFIX) \
	py_contour_style.$(OBJ_SUFFIX)

PY_CONTOUR_LEVELS_OBJS = \
	contour_levels.$(OBJ_SUFFIX) \
	py_contour_levels.$(OBJ_SUFFIX)

PY_PEAK_LIST_OBJS = \
	method.$(OBJ_SUFFIX) \
	peak.$(OBJ_SUFFIX) \
	peak_list.$(OBJ_SUFFIX) \
	symbol.$(OBJ_SUFFIX) \
	py_peak.$(OBJ_SUFFIX) \
	py_peak_list.$(OBJ_SUFFIX)

PY_WIN_PEAK_LIST_OBJS = \
	method.$(OBJ_SUFFIX) \
	peak.$(OBJ_SUFFIX) \
	peak_list.$(OBJ_SUFFIX) \
	symbol.$(OBJ_SUFFIX) \
	win_peak_list.$(OBJ_SUFFIX) \
	py_peak.$(OBJ_SUFFIX) \
	py_peak_list.$(OBJ_SUFFIX) \
	py_win_peak_list.$(OBJ_SUFFIX)

PY_PEAK_CLUSTER_OBJS = \
	method.$(OBJ_SUFFIX) \
	peak.$(OBJ_SUFFIX) \
	peak_cluster.$(OBJ_SUFFIX) \
	symbol.$(OBJ_SUFFIX) \
	py_peak.$(OBJ_SUFFIX) \
	py_peak_cluster.$(OBJ_SUFFIX)

all: ContourFile.$(DYLIB_SUFFIX) ContourStyle.$(DYLIB_SUFFIX) ContourLevels.$(DYLIB_SUFFIX) \
     PeakList.$(DYLIB_SUFFIX) SliceFile.$(DYLIB_SUFFIX) WinPeakList.$(DYLIB_SUFFIX) \
     PeakCluster.$(DYLIB_SUFFIX)

ContourFile.$(DYLIB_SUFFIX): global_contour_objects global_py_block_file_objects global_py_store_file_objects global_py_drawing_objects $(PY_CONTOUR_FILE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(ALL_LIB_FLAGS) \
		$(GLOBAL_CONTOUR_OBJS) $(GLOBAL_PY_BLOCK_FILE_OBJS) \
		$(GLOBAL_PY_STORE_FILE_OBJS) $(GLOBAL_PY_DRAWING_OBJS) \
		$(PY_CONTOUR_FILE_OBJS) $(ALL_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

SliceFile.$(DYLIB_SUFFIX): global_py_block_file_objects global_py_drawing_objects $(PY_SLICE_FILE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(ALL_LIB_FLAGS) \
		$(GLOBAL_PY_BLOCK_FILE_OBJS) $(GLOBAL_PY_DRAWING_OBJS) \
		$(PY_SLICE_FILE_OBJS) $(ALL_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

ContourStyle.$(DYLIB_SUFFIX): global_py_util_objects $(PY_CONTOUR_STYLE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_CONTOUR_STYLE_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

ContourLevels.$(DYLIB_SUFFIX): global_py_util_objects $(PY_CONTOUR_LEVELS_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_PY_UTIL_OBJS) $(PY_CONTOUR_LEVELS_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

PeakList.$(DYLIB_SUFFIX): global_nonlinear_model_objects global_py_block_file_objects $(PY_PEAK_LIST_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(GLOBAL_NONLINEAR_MODEL_OBJS) $(GLOBAL_PY_BLOCK_FILE_OBJS) \
		$(PY_PEAK_LIST_OBJS) $(MATH_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

WinPeakList.$(DYLIB_SUFFIX): global_nonlinear_model_objects global_py_block_file_objects global_py_drawing_objects $(PY_WIN_PEAK_LIST_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(ALL_LIB_FLAGS) \
		$(GLOBAL_NONLINEAR_MODEL_OBJS) $(GLOBAL_PY_BLOCK_FILE_OBJS) $(GLOBAL_PY_DRAWING_OBJS) \
		$(PY_WIN_PEAK_LIST_OBJS) $(ALL_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

PeakCluster.$(DYLIB_SUFFIX): global_nonlinear_model_objects global_py_block_file_objects global_py_drawing_objects $(PY_PEAK_CLUSTER_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(ALL_LIB_FLAGS) \
		$(GLOBAL_NONLINEAR_MODEL_OBJS) $(GLOBAL_PY_BLOCK_FILE_OBJS) $(GLOBAL_PY_DRAWING_OBJS) \
		$(PY_PEAK_CLUSTER_OBJS) $(ALL_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

py_contour_file.$(OBJ_SUFFIX): py_contour_file.c
	$(CC) $(CC_FLAGS) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_slice_file.$(OBJ_SUFFIX): py_slice_file.c
	$(CC) $(CC_FLAGS) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_contour_style.$(OBJ_SUFFIX): py_contour_style.c
	$(CC) $(CC_FLAGS) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_contour_levels.$(OBJ_SUFFIX): py_contour_levels.c
	$(CC) $(CC_FLAGS) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_peak_list.$(OBJ_SUFFIX): py_peak_list.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_peak.$(OBJ_SUFFIX): py_peak.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_win_peak_list.$(OBJ_SUFFIX): py_win_peak_list.c
	$(CC) $(CC_FLAGS) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_peak_cluster.$(OBJ_SUFFIX): py_peak_cluster.c
	$(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

.c.$(OBJ_SUFFIX):;  $(CC) $(CC_FLAGS) $(GLOBAL_INCLUDE_FLAGS) $(CO_NAME)
 
global_py_util_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_util_objects

global_py_mem_cache_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_mem_cache_objects

global_py_block_file_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_block_file_objects

global_py_store_file_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_store_file_objects

global_py_shape_file_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_shape_file_objects

global_contour_objects:
	cd $(GLOBAL_DIR) && $(MAKE) contour_objects

global_drawing_objects:
	cd $(GLOBAL_DIR) && $(MAKE) drawing_objects

global_py_drawing_objects:
	cd $(GLOBAL_DIR) && $(MAKE) py_drawing_objects

global_nonlinear_model_objects:
	cd $(GLOBAL_DIR) && $(MAKE) nonlinear_model_objects

clean:
	$(RM) *.$(OBJ_SUFFIX)

realclean: clean
	$(RM) *.$(DYLIB_SUFFIX)

