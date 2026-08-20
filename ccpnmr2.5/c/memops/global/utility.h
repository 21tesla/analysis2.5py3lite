
/*
======================COPYRIGHT/LICENSE START==========================

utility.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_utility
#define _incl_utility

#include "macros.h"
#include "types.h"

#define  IS_BIG_ENDIAN  0
#define  IS_LITTLE_ENDIAN  1

/* below should be used only on Azara-written data */
extern Bool determine_swapped
	(void);

extern CcpnStatus endian_fwrite
	(char *x, int n, FILE *fp);

extern Bool is_big_endian
	(void);

extern Bool is_little_endian
	(void);

extern int endianess
	(void);

extern Bool file_exists
	(CcpnString file);

extern Bool equal_strings
	(CcpnString string1, CcpnString string2);

extern CcpnStatus copy_string
	(CcpnString *string1, CcpnString string2, CcpnString error_msg);

extern Bool empty_string
	(CcpnString string);

extern void print_integer_array
	(int n, int *array);

extern void swap_bytes
	(char *ptr, int nbytes);

extern void float_words
	(float *ptr, int nwords);

extern void int_words
	(int *ptr, int nwords);

extern int log_2
	(int n);

extern int floor_power_of_2
	(int n);

extern int ceil_power_of_2
	(int n);

extern int greatest_common_factor
	(int m, int n);

extern Bool no_abort_func
	(void);

extern YesNo yes_func
	(CcpnString msg, CcpnString yes, CcpnString no, CcpnString cancel);

extern YesNo no_func
	(CcpnString msg, CcpnString yes, CcpnString no, CcpnString cancel);

extern YesNo cancel_func
	(CcpnString msg, CcpnString yes, CcpnString no, CcpnString cancel);

extern CcpnStatus start_screen_print
	(CcpnString error_msg);

extern void print_screen_message
	(CcpnString message);

extern void end_screen_print
	(void);

extern void start_no_timer
	(CcpnString message);

extern int update_no_timer
	(float fraction);

extern void stop_no_timer
	(CcpnStatus status);

extern void start_screen_timer
	(CcpnString message);

extern int update_screen_timer
	(float fraction);

extern void stop_screen_timer
	(CcpnStatus status);

extern CcpnStatus get_integers
	(int n, int *array, char *string);

extern CcpnStatus get_floats
	(int n, float *array, char *string);

extern int get_max_integers
	(int nmax, int *array, char *string);

extern int get_max_floats
	(int nmax, float *array, char *string);

extern CcpnStatus get_some_integers
	(int n, int m, int *skip, int *array, char *string);

extern CcpnStatus get_some_floats
	(int n, int m, int *skip, float *array, char *string);

extern int fseek_absolute_bytes(FILE *fp, Integer nbytes);

#endif /* _incl_utility */
