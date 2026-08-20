
/*
======================COPYRIGHT/LICENSE START==========================

types.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_types
#define _incl_types

#include "macros.h"

#ifdef Bool
#undef Bool
#endif

typedef char  *CcpnString;

typedef enum { CCPN_OK, CCPN_ERROR, OTHER } CcpnStatus;
typedef enum { NEITHER, YES, NO, CANCEL } YesNo;

typedef void Generic; /* TBD: change */
typedef Generic *Generic_ptr; /* TBD: change */
typedef void *CcpnGenericPtr;
typedef int Bool; /* TBD: change all Bool to CcpnBool */
typedef int CcpnBool;
typedef char Line[LINE_SIZE+1];
typedef char Long_line[LONG_LINE_SIZE+1];
typedef char Param_value[LINE_SIZE+1];

typedef long long int Integer;
typedef float Float;

typedef void (*Transform)(float *, float *);

typedef int (*Compare_func)(Generic_ptr, Generic_ptr);
typedef Bool (*Equal_func)(Generic_ptr, Generic_ptr);
typedef CcpnStatus (*Check_func)(Generic_ptr);
typedef void (*Data_func)(Generic_ptr);
typedef Bool (*Truth_func)(Generic_ptr);
typedef int (*Search_func)(int);
typedef unsigned long Hash_value;
typedef Hash_value (*Hash_func)(Generic_ptr);
typedef Bool (*Abort_func)(void);

typedef YesNo (*YesNo_func)(CcpnString msg, CcpnString yes, CcpnString no, CcpnString cancel);

typedef struct
{
    CcpnStatus (*start_print)(CcpnString error_msg);
    void (*print_message)(CcpnString message);
    void (*end_print)(void);
}   Print_funcs;

typedef struct
{
    void (*start_timer)(CcpnString message);
    int (*update_timer)(float fraction);
    void (*stop_timer)(CcpnStatus status);
}   Timer_funcs;

#define  SCREEN_DISPLAY		0
#define  PAPER_DISPLAY		1

typedef struct
{
    int display_medium;
    Generic_ptr data;
    void (*start_draw)(Generic_ptr data);
    void (*end_draw)(void);
    void (*new_draw_range)(float x0, float y0, float x1, float y1, Bool clip);
    void (*draw_line)(float x0, float y0, float x1, float y1);
    void (*draw_text)(CcpnString text, float x, float y, float a, float b);
    void (*set_draw_color)(int color);
    void (*set_draw_font)(CcpnString name, int size);
    void (*set_line_style)(int line_style);
}   Draw_funcs;

#endif /* _incl_types */
