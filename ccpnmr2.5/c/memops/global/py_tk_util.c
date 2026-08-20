
/*
======================COPYRIGHT/LICENSE START==========================

py_tk_util.c: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
*/

#include "py_tk_util.h"

/*
   The Python Tkinter interface to Tk doesn't provide public access
   to its C types.  This typedef is from _tkinter.c.
*/

//struct TkappObject
//{
//    PyObject_HEAD
//    Tcl_Interp *interp;
//};

typedef struct {
    PyObject_HEAD
    Tcl_Interp *interp;
    int wantobjects;
    int threaded; /* True if tcl_platform[threaded] */
    Tcl_ThreadId thread_id;
    int dispatching;
    /* We cannot include tclInt.h, as this is internal.
       So we cache interesting types here. */
    const Tcl_ObjType *OldBooleanType;
    const Tcl_ObjType *BooleanType;
    const Tcl_ObjType *ByteArrayType;
    const Tcl_ObjType *DoubleType;
    const Tcl_ObjType *IntType;
    const Tcl_ObjType *WideIntType;
    const Tcl_ObjType *BignumType;
    const Tcl_ObjType *ListType;
    const Tcl_ObjType *ProcBodyType;
    const Tcl_ObjType *StringType;
} TkappObject;

#define Tkapp_Interp(v) (((TkappObject *) (v))->interp)

Tcl_Interp *get_tcl_interp(PyObject *widget, CcpnString error_msg)
{
    PyObject *tkapp = NULL;
    Tcl_Interp *tcl_interp = NULL;

    if (!PyObject_HasAttrString(widget, "tk"))
    {
        sprintf(error_msg, "argument not a Tk widget");
    }
    else
    {
	tkapp = PyObject_GetAttrString(widget, "tk");

	if (tkapp == Py_None)
        sprintf(error_msg, "widget tk is None");
	else
//	    tcl_interp = ((struct TkappObject *) tkapp)->interp;
	    tcl_interp = Tkapp_Interp(tkapp);
    }

    Py_XDECREF(tkapp);

    return tcl_interp;
}

Tk_Window get_tk_window(PyObject *widget, Tcl_Interp *tcl_interp,
						CcpnString error_msg)
{
    PyObject *o = NULL;
    char *path;
    Tk_Window tk_main_win, tk_display_win = NULL;
    int tk_count = -1;

    if (!PyObject_HasAttrString(widget, "_w"))
    {
        sprintf(error_msg, "argument not a Tk widget");
    }
    else 
    {
        o = PyObject_GetAttrString(widget, "_w");
        if (o == Py_None)
	    {
            sprintf(error_msg, "widget _w is None");
	    }
        else
        {
    //	    Tcl_Interp *TCLinterpLocal = NULL;
    //        TCLinterpLocal = Tcl_CreateInterp();
    //	    tk_main_win = (Tk_Window) Tk_MainWindow(TCLinterpLocal);

            tk_main_win = (Tk_Window) Tk_MainWindow(tcl_interp);
            tk_count = (int) Tk_GetNumMainWindows();
            if (tk_main_win)
            {
                path = PyString_AsString(o);
                tk_display_win = Tk_NameToWindow(tcl_interp, path, tk_main_win);
            }
            else
            {
//                path = PyString_AsString(o);
                sprintf(error_msg, "could not get main window, possibly different Python Tcl/Tk and Analysis Tcl/Tk (%d windows).", tk_count);
            }
        }
    }

    Py_XDECREF(o);

    return tk_display_win;
}

