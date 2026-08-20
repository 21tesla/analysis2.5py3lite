include ../../environment.txt
# above gives CC, LINK, MAKE, CC_FLAGS, LINKFLAGS, ..._INCLUDE_FLAGS, ..._LIB_FLAGS, ..._LIB

UTIL_OBJS = \
	utility.$(OBJ_SUFFIX)

PY_UTIL_OBJS = \
	$(UTIL_OBJS) \
	python_util.$(OBJ_SUFFIX)

COLOR_OBJS = \
	color.$(OBJ_SUFFIX)

SORTS_OBJS = \
	sorts.$(OBJ_SUFFIX)

LINALG_OBJS = \
	eigenvalue.$(OBJ_SUFFIX) \
	linalg.$(OBJ_SUFFIX)

GEOMETRY_OBJS = \
	geometry.$(OBJ_SUFFIX)

CLIPPING_OBJS = \
	clipping.$(OBJ_SUFFIX)

PY_MEM_CACHE_OBJS = \
	$(PY_UTIL_OBJS) \
	hash_list.$(OBJ_SUFFIX) \
	hash_table.$(OBJ_SUFFIX) \
	mem_cache.$(OBJ_SUFFIX) \
	mutex.$(OBJ_SUFFIX) \
	py_mem_cache.$(OBJ_SUFFIX)

PY_BLOCK_FILE_OBJS = \
	$(PY_MEM_CACHE_OBJS) \
	block_file.$(OBJ_SUFFIX) \
	shape_file.$(OBJ_SUFFIX) \
	int_array.$(OBJ_SUFFIX) \
	py_block_file.$(OBJ_SUFFIX) \
	py_shape_file.$(OBJ_SUFFIX)

CONTOUR_OBJS = \
	contourer.$(OBJ_SUFFIX) \
	list.$(OBJ_SUFFIX)

PY_GL_HANDLER_OBJS = \
	$(PY_UTIL_OBJS) \
	$(CLIPPING_OBJS) \
	gl_handler.$(OBJ_SUFFIX) \
	py_gl_handler.$(OBJ_SUFFIX) \
	py_tk_util.$(OBJ_SUFFIX) \
	w32func.$(OBJ_SUFFIX)

PY_TK_HANDLER_OBJS = \
	$(PY_UTIL_OBJS) \
	$(CLIPPING_OBJS) \
	$(COLOR_OBJS) \
	tk_handler.$(OBJ_SUFFIX) \
	py_tk_handler.$(OBJ_SUFFIX) \
	py_tk_util.$(OBJ_SUFFIX) \
	w32func.$(OBJ_SUFFIX)

PY_PDF_HANDLER_OBJS = \
	$(PY_UTIL_OBJS) \
	$(CLIPPING_OBJS) \
	pdf_handler.$(OBJ_SUFFIX) \
	py_pdf_handler.$(OBJ_SUFFIX)

PY_PS_HANDLER_OBJS = \
	$(PY_UTIL_OBJS) \
	$(CLIPPING_OBJS) \
	ps_handler.$(OBJ_SUFFIX) \
	py_ps_handler.$(OBJ_SUFFIX)

PY_STORE_HANDLER_OBJS = \
	$(PY_UTIL_OBJS) \
	store_handler.$(OBJ_SUFFIX) \
	py_store_handler.$(OBJ_SUFFIX)

PY_STORE_FILE_OBJS = \
	store_file.$(OBJ_SUFFIX) \
	py_store_file.$(OBJ_SUFFIX)

PY_SHAPE_FILE_OBJS = \
	$(PY_UTIL_OBJS) \
	shape_file.$(OBJ_SUFFIX) \
	hash_list.$(OBJ_SUFFIX) \
	hash_table.$(OBJ_SUFFIX) \
	int_array.$(OBJ_SUFFIX) \
	py_shape_file.$(OBJ_SUFFIX)

DRAWING_OBJS = \
	$(COLOR_OBJS) \
	$(CLIPPING_OBJS) \
	gl_handler.$(OBJ_SUFFIX) \
	pdf_handler.$(OBJ_SUFFIX) \
	ps_handler.$(OBJ_SUFFIX) \
	store_handler.$(OBJ_SUFFIX) \
	tk_handler.$(OBJ_SUFFIX) \
	w32func.$(OBJ_SUFFIX)

PY_DRAWING_OBJS = \
	$(DRAWING_OBJS) \
	py_draw_handler.$(OBJ_SUFFIX) \
	py_gl_handler.$(OBJ_SUFFIX) \
	py_pdf_handler.$(OBJ_SUFFIX) \
	py_ps_handler.$(OBJ_SUFFIX) \
	py_store_handler.$(OBJ_SUFFIX) \
	py_tk_handler.$(OBJ_SUFFIX) \
	py_tk_util.$(OBJ_SUFFIX)

RANDOM_OBJS = \
	random.$(OBJ_SUFFIX)

DIAG_OBJS = \
	diag_dbl.$(OBJ_SUFFIX)

NONLINEAR_MODEL_OBJS = \
	nonlinear_model.$(OBJ_SUFFIX) \
	gauss_jordan.$(OBJ_SUFFIX)

FIT_OBJS = \
	$(NONLINEAR_MODEL_OBJS) \
	fit.$(OBJ_SUFFIX) \
	cpmg.$(OBJ_SUFFIX) \
	fit1d.$(OBJ_SUFFIX) \
	gamma.$(OBJ_SUFFIX) \
	line_fit.$(OBJ_SUFFIX) \
	random.$(OBJ_SUFFIX)

PY_FIT_OBJS = \
	$(FIT_OBJS) \
	$(PY_UTIL_OBJS) \
	py_fit.$(OBJ_SUFFIX)

TCLTK_INCLUDE_FLAGS = $(PYTHON_INCLUDE_FLAGS) $(X11_INCLUDE_FLAGS) \
	$(TCL_INCLUDE_FLAGS) $(TK_INCLUDE_FLAGS)

TCLTK_LIB_FLAGS = $(X11_LIB_FLAGS) $(TCL_LIB_FLAGS) $(TK_LIB_FLAGS)

TCLTK_LIB = $(TK_LIB) $(X11_LIB) $(MATH_LIB)

ALL_INCLUDE_FLAGS = $(PYTHON_INCLUDE_FLAGS) $(X11_INCLUDE_FLAGS) \
	$(TCL_INCLUDE_FLAGS) $(TK_INCLUDE_FLAGS) $(GL_INCLUDE_FLAGS)

ALL_LIB_FLAGS = $(X11_LIB_FLAGS) $(TCL_LIB_FLAGS) $(TK_LIB_FLAGS) $(GL_LIB_FLAGS)

ALL_LIB = $(GL_LIB) $(TK_LIB) $(X11_LIB) $(MATH_LIB)

all: BlockFile.$(DYLIB_SUFFIX) GlHandler.$(DYLIB_SUFFIX) MemCache.$(DYLIB_SUFFIX) PdfHandler.$(DYLIB_SUFFIX) PsHandler.$(DYLIB_SUFFIX) FitMethod.$(DYLIB_SUFFIX) TkHandler.$(DYLIB_SUFFIX) StoreFile.$(DYLIB_SUFFIX) StoreHandler.$(DYLIB_SUFFIX) ShapeFile.$(DYLIB_SUFFIX)

util_objects: $(UTIL_OBJS)
py_util_objects: $(PY_UTIL_OBJS)
color_objects: $(COLOR_OBJS)
sorts_objects: $(SORTS_OBJS)
linalg_objects: $(LINALG_OBJS)
geometry_objects: $(GEOMETRY_OBJS)
clipping_objects: $(CLIPPING_OBJS)
py_mem_cache_objects: $(PY_MEM_CACHE_OBJS)
py_block_file_objects: $(PY_BLOCK_FILE_OBJS)
py_store_file_objects: $(PY_STORE_FILE_OBJS)
py_shape_file_objects: $(PY_SHAPE_FILE_OBJS)
contour_objects: $(CONTOUR_OBJS)
drawing_objects: $(DRAWING_OBJS)
py_drawing_objects: $(PY_DRAWING_OBJS)
random_objects: $(RANDOM_OBJS)
diag_objects: $(DIAG_OBJS)
nonlinear_model_objects: $(NONLINEAR_MODEL_OBJS)
py_fit_objects: $(PY_FIT_OBJS)

BlockFile.$(DYLIB_SUFFIX): $(PY_BLOCK_FILE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(PY_BLOCK_FILE_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

GlHandler.$(DYLIB_SUFFIX): $(PY_GL_HANDLER_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(ALL_LIB_FLAGS) \
		$(PY_GL_HANDLER_OBJS) $(ALL_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

TkHandler.$(DYLIB_SUFFIX): $(PY_TK_HANDLER_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(TK_LIB_FLAGS) $(X11_LIB_FLAGS) \
		$(PY_TK_HANDLER_OBJS) $(TK_LIB) $(X11_LIB) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

PdfHandler.$(DYLIB_SUFFIX): $(PY_PDF_HANDLER_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(PY_PDF_HANDLER_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

PsHandler.$(DYLIB_SUFFIX): $(PY_PS_HANDLER_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(PY_PS_HANDLER_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

StoreFile.$(DYLIB_SUFFIX): $(PY_UTIL_OBJS) $(PY_STORE_FILE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(PY_UTIL_OBJS) $(PY_STORE_FILE_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

StoreHandler.$(DYLIB_SUFFIX): $(PY_STORE_HANDLER_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(PY_STORE_HANDLER_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

ShapeFile.$(DYLIB_SUFFIX): $(PY_SHAPE_FILE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) \
		$(PY_SHAPE_FILE_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

MemCache.$(DYLIB_SUFFIX): $(PY_MEM_CACHE_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(PY_MEM_CACHE_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

FitMethod.$(DYLIB_SUFFIX): $(PY_FIT_OBJS)
	$(LINK) $(LINK_FLAGS) $(SHARED_FLAGS) $(OUT_NAME) $(PY_FIT_OBJS) $(PYTHON_LIB_FLAGS) $(PYTHON_LIB)

py_mem_cache.$(OBJ_SUFFIX): py_mem_cache.c 
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_block_file.$(OBJ_SUFFIX): py_block_file.c 
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_store_file.$(OBJ_SUFFIX): py_store_file.c 
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_shape_file.$(OBJ_SUFFIX): py_shape_file.c 
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_fit.$(OBJ_SUFFIX): py_fit.c 
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

python_util.$(OBJ_SUFFIX): python_util.c 
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_draw_handler.$(OBJ_SUFFIX): py_draw_handler.c
	$(CC) $(CC_FLAGS) $(IGNORE_GL_FLAG) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_gl_handler.$(OBJ_SUFFIX): py_gl_handler.c
	$(CC) $(CC_FLAGS) $(IGNORE_GL_FLAG) $(GL_FLAG) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

gl_handler.$(OBJ_SUFFIX): gl_handler.c
	$(CC) $(CC_FLAGS) $(IGNORE_GL_FLAG) $(GLUT_FLAG) $(ALL_INCLUDE_FLAGS) $(CO_NAME)

py_tk_handler.$(OBJ_SUFFIX): py_tk_handler.c
	$(CC) $(CC_FLAGS) $(TCLTK_INCLUDE_FLAGS) $(CO_NAME)

tk_handler.$(OBJ_SUFFIX): tk_handler.c
	$(CC) $(CC_FLAGS) $(TCLTK_INCLUDE_FLAGS) $(CO_NAME)

py_tk_util.$(OBJ_SUFFIX): py_tk_util.c
	$(CC) $(CC_FLAGS) $(TCLTK_INCLUDE_FLAGS) $(CO_NAME)

py_ps_handler.$(OBJ_SUFFIX): py_ps_handler.c
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_store_handler.$(OBJ_SUFFIX): py_store_handler.c
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

py_pdf_handler.$(OBJ_SUFFIX): py_pdf_handler.c
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDE_FLAGS) $(CO_NAME)

.c.$(OBJ_SUFFIX):;  $(CC) $(CC_FLAGS) $(CO_NAME)
 
clean:
	$(RM) *.$(OBJ_SUFFIX)
 
realclean: clean
	$(RM) *.$(DYLIB_SUFFIX)

