
/*
======================COPYRIGHT/LICENSE START==========================

utility.c: Part of the CcpNmr Analysis program

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
#include "utility.h"

Bool determine_swapped()
{
    Bool swapped = CCPN_FALSE;

#ifdef BIG_ENDIAN_DATA
    if (is_little_endian())
	swapped = CCPN_TRUE;
#endif /* BIG_ENDIAN_DATA */

#ifdef LITTLE_ENDIAN_DATA
    if (is_big_endian())
	swapped = CCPN_TRUE;
#endif /* LITTLE_ENDIAN_DATA */

    return  swapped;
}

CcpnStatus endian_fwrite(char *x, int n, FILE *fp)
{
    CcpnStatus status = CCPN_OK;
    Bool swapped = determine_swapped();

    if (swapped)
	swap_bytes(x, BYTES_PER_WORD*n);

    if (FWRITE(x, BYTES_PER_WORD, n, fp))
	status = CCPN_ERROR;

    if (swapped)
	swap_bytes(x, BYTES_PER_WORD*n);

    return  status;
}

/* below only handles four byte case, not sure about rest */
Bool is_big_endian()
{
    static char c[4] = { 0, 0, 0, 1 };
    int x;

    x = *((int *) &c);

    if (x == 1)
	return  CCPN_TRUE;
    else
	return  CCPN_FALSE;
}

Bool is_little_endian()
{
    return !is_big_endian();
}

int endianess()
{
    if (is_big_endian())
	return  IS_BIG_ENDIAN;
    else
	return  IS_LITTLE_ENDIAN;
}

Bool file_exists(CcpnString file)
{
    FILE *fp;

    if (OPEN_FOR_READING(fp, file))
    {
	return  CCPN_FALSE;
    }
    else
    {
	fclose(fp);

	return  CCPN_TRUE;
    }
}

Bool equal_strings(CcpnString string1, CcpnString string2)
{
    if (!string1 && string2)
	return  CCPN_FALSE;

    if (!string2 && string1)
	return  CCPN_FALSE;

    if (strcmp(string1, string2))
	return  CCPN_FALSE;
    else
	return  CCPN_TRUE;
}

CcpnStatus copy_string(CcpnString *string1, CcpnString string2, CcpnString error_msg)
{
    int n = strlen(string2) + 1;
    char *s;

    sprintf(error_msg, "allocating memory for string copy");
    MALLOC(s, char, n);

    strcpy(s, string2);
    *string1 = s;

    return  CCPN_OK;
}

Bool empty_string(CcpnString string)
{
    int i, n = strlen(string);

    for (i = 0; i < n; i++)
    {
	if (!isspace(string[i]))
	    return  CCPN_FALSE;
    }

    return  CCPN_TRUE;
}

/* print_integer_array assumes array components should have +1 added */

void print_integer_array(int n, int *array)
{
    int i;

    printf("(");

    for (i = 0; i < n; i++)
    {
        if (i > 0)
            printf(", ");

        printf("%d", array[i]+1);
    }

    printf(")");
}

/* swap_bytes assumes that there are 4 bytes per word */

void swap_bytes(char *ptr, int nbytes)
{
    char *ptr2;

    for (ptr2 = ptr+nbytes; ptr < ptr2; ptr += 4)
    {
	SWAP(ptr[0], ptr[3], char);
	SWAP(ptr[1], ptr[2], char);
    }
}

/* float_words takes integer data and converts it to floating point */

void float_words(float *ptr, int nwords)
{
    float *ptr2;

    for (ptr2 = ptr+nwords; ptr < ptr2; ptr++)
	*ptr = *((int *) ptr);
}

/* int_words takes floating point data and converts it to integer by rounding */

void int_words(int *ptr, int nwords)
{
    int *ptr2;
    float x;

    for (ptr2 = ptr+nwords; ptr < ptr2; ptr++)
    {
	x = *((float *) ptr);
	*ptr = NEAREST_INTEGER(x);
    }
}

int log_2(int n)  /* assumes that n is a power of 2 */
{
    int m = -1;

    while (n > 0)
    {
	n >>= 1;  m++;
    }

    return  m;
}

int floor_power_of_2(int n)
{
    int m = -1;

    if (n < 1)  return  0;  /* arbitrary */

    while (n > 0)
    {
	n >>= 1;  m++;
    }

    return  (1 << m);
}

int ceil_power_of_2(int n)
{
    int m = -1, p, q = n;

    if (q < 1)  return  0;  /* arbitrary */

    while (q > 0)
    {
	q >>= 1;  m++;
    }

    p = 1 << m;

    if (p < n)
	return  2*p;
    else
	return  p;
}

int greatest_common_factor(int m, int n)
{
    int p;

    if (m == 0)
	return  n;

    if (n == 0)
	return  m;

    m = ABS(m);
    n = ABS(n);

    if (m > n)
	SWAP(m, n, int);

    p = n % m;
    while (p != 0)
    {
	n = m;
	m = p;
	p = n % m;
    }

    return  m;
}

Bool no_abort_func()
{
    return  CCPN_FALSE;
}

YesNo yes_func(CcpnString msg, CcpnString yes, CcpnString no, CcpnString cancel)
{
    return  YES;
}

YesNo no_func(CcpnString msg, CcpnString yes, CcpnString no, CcpnString cancel)
{
    return  NO;
}

YesNo cancel_func(CcpnString msg, CcpnString yes, CcpnString no, CcpnString cancel)
{
    return  CANCEL;
}

CcpnStatus start_screen_print(CcpnString error_msg)
{
    return  CCPN_OK;
}

void print_screen_message(CcpnString message)
{
    printf("%s", message);
}

void end_screen_print()
{
}

void start_no_timer(CcpnString message)
{
}

int update_no_timer(float fraction)
{
    return  CONTINUE;
}

void stop_no_timer(CcpnStatus status)
{
}

void start_screen_timer(CcpnString message)
{
    printf("%s ...\n", message);
}

int update_screen_timer(float fraction)
{
    int percent = 100.0 * fraction;
    printf("\t%3d%% done\n", percent);

    return  CONTINUE;
}

void stop_screen_timer(CcpnStatus status)
{
    if (status == CCPN_OK)
	printf("... Finished\n");
}

#define  SEPARATOR  " \t"

CcpnStatus get_integers(int n, int *array, char *string)
{
    int i;
    char *ptr1, *ptr2;
    Line copy;

    strncpy(copy, string, LINE_SIZE);
    ptr1 = copy;
    for (i = 0; i < n; i++)
    {
	if (!(ptr2 = strtok(ptr1, SEPARATOR)) ||
		(sscanf(ptr2, "%d", array+i) != 1))
	    return  CCPN_ERROR;

	ptr1 = NULL;
    }

    return  CCPN_OK;
}

CcpnStatus get_floats(int n, float *array, char *string)
{
    int i;
    char *ptr1, *ptr2;
    Line copy;

    strncpy(copy, string, LINE_SIZE);
    ptr1 = copy;
    for (i = 0; i < n; i++)
    {
	if (!(ptr2 = strtok(ptr1, SEPARATOR)) ||
		(sscanf(ptr2, "%f", array+i) != 1))
	    return  CCPN_ERROR;

	ptr1 = NULL;
    }

    return  CCPN_OK;
}

int get_max_integers(int nmax, int *array, char *string)
{
    int n;
    char *ptr1, *ptr2;
    Line copy;

    strncpy(copy, string, LINE_SIZE);
    ptr1 = copy;
    for (n = 0; n < nmax; n++)
    {
	if (!(ptr2 = strtok(ptr1, SEPARATOR)) ||
		(sscanf(ptr2, "%d", array+n) != 1))
	    return  n;

	ptr1 = NULL;
    }

    return  n;
}

int get_max_floats(int nmax, float *array, char *string)
{
    int n;
    char *ptr1, *ptr2;
    Line copy;

    strncpy(copy, string, LINE_SIZE);
    ptr1 = copy;
    for (n = 0; n < nmax; n++)
    {
	if (!(ptr2 = strtok(ptr1, SEPARATOR)) ||
		(sscanf(ptr2, "%f", array+n) != 1))
	    return  n;

	ptr1 = NULL;
    }

    return  n;
}

static Bool skip_this(int ind, int n, int *skip)
{
    int i;

    for (i = 0; i < n; i++)
	if (skip[i] == ind)
	    return  CCPN_TRUE;

    return  CCPN_FALSE;
}

CcpnStatus get_some_integers(int n, int m, int *skip, int *array, char *string)
{
    int i;
    char *ptr1, *ptr2;
    Line copy;

    strncpy(copy, string, LINE_SIZE);
    ptr1 = copy;
    for (i = 0; i < n; i++)
    {
	if (!(ptr2 = strtok(ptr1, SEPARATOR)))
	    return  CCPN_ERROR;

	if (skip_this(i, m, skip))
	    array[i] = 0;
	else if (sscanf(ptr2, "%d", array+i) != 1)
	    return  CCPN_ERROR;

	ptr1 = NULL;
    }

    return  CCPN_OK;
}

CcpnStatus get_some_floats(int n, int m, int *skip, float *array, char *string)
{
    int i;
    char *ptr1, *ptr2;
    Line copy;

    strncpy(copy, string, LINE_SIZE);
    ptr1 = copy;
    for (i = 0; i < n; i++)
    {
	if (!(ptr2 = strtok(ptr1, SEPARATOR)))
	    return  CCPN_ERROR;

	if (skip_this(i, m, skip))
	    array[i] = 0;
	else if (sscanf(ptr2, "%f", array+i) != 1)
	    return  CCPN_ERROR;

	ptr1 = NULL;
    }

    return  CCPN_OK;
}

#define  PADDING  "xyabcdzw"

static char padding[] = PADDING;
static int npadding = 8;

void *ccpn_malloc(size_t size)
{
    char *x;
    x = (char *) malloc(size+npadding);

    strncpy(x, padding, npadding);

    return (void *) (x+npadding);
}

void *ccpn_realloc(void *ptr, size_t size)
{
    char *x = ((char *) ptr) - npadding, *y;
    y = (char *) realloc(x, size+npadding);

    return (void *) (y+npadding);
}

void ccpn_free(void *ptr)
{
    char *x = ((char *) ptr) - npadding;

    if (strncmp(x, padding, npadding))
	printf("TRAMPLE: %p\n", ptr);
}

int fseek_absolute_bytes(FILE *fp, Integer nbytes)
{
    int result;
    Bool looped = CCPN_FALSE;
    long int nb, BIGN = 1<<30;

    while (nbytes > BIGN)
    {
	if (looped)
            result = FSEEK_RELATIVE_BYTES(fp, BIGN);
	else
            result = FSEEK_ABSOLUTE_BYTES(fp, BIGN);

        if (result)
	    return result;

        nbytes -= BIGN;
	looped = CCPN_TRUE;
    }

    nb = nbytes;

    if (looped)
        return FSEEK_RELATIVE_BYTES(fp, nb);
    else
        return FSEEK_ABSOLUTE_BYTES(fp, nb);
}

