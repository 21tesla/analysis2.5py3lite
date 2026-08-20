
/*
======================COPYRIGHT/LICENSE START==========================

macros.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_macros
#define _incl_macros

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#ifndef WIN64
#include <sys/types.h>
#endif /* WIN32 */

#ifdef WIN64
#define  DIRECTORY_SYMBOL  '\\'
extern int truncate(char *file, long length);
extern void sleep(int seconds);
#else /* WIN32 */
#define  DIRECTORY_SYMBOL  '/'
#endif /* WIN32 */

extern void print_error_message( void );

/* normal mode */
#define  DEBUG_BLOCK(x)  
/* debug mode */
/*
#define  DEBUG_BLOCK(x)  x
*/

#define  DEBUG_FILE  "log.txt"

#define  DEBUG_CODE_FMT(fmt, val)  DEBUG_BLOCK({  FILE *fp_dbg = fopen(DEBUG_FILE, "a"); fprintf(fp_dbg, "%s %d " fmt "\n", __FILE__, __LINE__, val); fclose(fp_dbg);   })

#define  DEBUG_CODE_FMT2(fmt, val1, val2)  DEBUG_BLOCK({  FILE *fp_dbg = fopen(DEBUG_FILE, "a"); fprintf(fp_dbg, "%s %d " fmt "\n", __FILE__, __LINE__, val1, val2); fclose(fp_dbg);   })

#define  DEBUG_CODE_HEX(msg, val)  DEBUG_CODE_FMT(msg " %x", val);

#define  DEBUG_CODE_HEX_INT(msg, val1, val2)  DEBUG_CODE_FMT2(msg " %x %d", val1, val2);

extern void *ccpn_malloc(size_t size);

extern void *ccpn_realloc(void *ptr, size_t size);

extern void ccpn_free(void *ptr);

#define  BYTES_PER_WORD  4

#define  CONTINUE  0
#define  ABORT  1

#define  BEGIN		0
#define  END		1

#define  ARGS(x)  ()

#define  LARGE_FLOAT  1.0e20

#ifndef  LINE_SIZE
#define  LINE_SIZE  120
#endif

#ifndef  LONG_LINE_SIZE
#define  LONG_LINE_SIZE  500
#endif

#ifndef  PI
#define  PI  3.14159265358979
#endif /* PI */

#define  TWOPI  (2.0 * PI)
#define  RADIAN_SCALE  (PI / 180.0)

#ifndef  HALF
#define  HALF  0.5
#endif /* HALF */

#define  SQRT2  (sqrt((double) 2.0))
#define  SQRT2_INV  (1.0 / SQRT2)
#define  LN2  (log(2.0))

#define  READ		"r"
#define  WRITE		"w"
#define  APPEND		"a"
#define  READ_UPDATE	"r+"
#define  WRITE_UPDATE	"w+"
#define  BINARY_READ	"rb"
#define  BINARY_WRITE	"wb"
#define  BINARY_READ_UPDATE	"r+b"
#define  BINARY_WRITE_UPDATE	"w+b"

/*
#define  R  'r'
#define  W  'w'
#define  A  'a'
*/

#define  END_OF_STRING  '\0'
#define  NEW_LINE  '\n'

#ifndef  CCPN_TRUE
#define  CCPN_TRUE   1
#endif

#ifndef  CCPN_FALSE
#define  CCPN_FALSE  0
#endif

#ifndef  SEEK_SET
#define  SEEK_SET  0
#define  SEEK_CUR  1
#define  SEEK_END  2
#endif

#ifndef  MIN
#define  MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef  MAX
#define  MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

#define  ABS(a)  (((a) < 0) ? -(a) : (a))

#define  SIGN(a)  (((a) < 0) ? -1 : 1)

#define  FLOOR(x) \
(((x) == (int) (x)) ? (int) (x) : ((x) >= 0  ?  (int) (x)  :  (- 1 - ((int) (-(x))))))

#define  CEILING(x)  (- FLOOR(-(x)))

#define  NEAREST_INTEGER(a)  (FLOOR((a) + 0.5))

#define  SWAP(a, b, type) \
	 {   type  I;  I = a;  a = b;  b = I;  }

#define  SWAP3(a, b, c, type) \
	 {   type  I;  I = a;  a = b;  b = c;  c = I;  }

#define  FLUSH  {   fflush(stdout);  fflush(stderr);   }

#define  OPEN_FOR_READING(fp, name) \
	 (((fp) = fopen((name), READ)) == NULL)

#define  OPEN_FOR_WRITING(fp, name) \
	 (((fp) = fopen((name), WRITE)) == NULL)

#define  OPEN_FOR_APPENDING(fp, name) \
	 (((fp) = fopen((name), APPEND)) == NULL)

#define  OPEN_FOR_BINARY_READING(fp, name) \
	 (((fp) = fopen((name), BINARY_READ)) == NULL)

#define  OPEN_FOR_BINARY_WRITING(fp, name) \
	 (((fp) = fopen((name), BINARY_WRITE)) == NULL)

#define  OPEN_FOR_BINARY_READ_UPDATE(fp, name) \
	 (((fp) = fopen((name), BINARY_READ_UPDATE)) == NULL)

#define  OPEN_FOR_BINARY_WRITE_UPDATE(fp, name) \
	 (((fp) = fopen((name), BINARY_WRITE_UPDATE)) == NULL)

#define  CHECK_OPEN_FOR_READING(fp, name) \
	 {   if (OPEN_FOR_READING(fp, name)) \
	     {   sprintf(error_msg, "opening '%s' for reading", name); \
		 return  CCPN_ERROR;   }   }

#define  CHECK_OPEN_FOR_BINARY_READING(fp, name) \
	 {   if (OPEN_FOR_BINARY_READING(fp, name)) \
	     {   sprintf(error_msg, "opening '%s' for reading", name); \
		 return  CCPN_ERROR;   }   }

#define  CHECK_OPEN_FOR_BINARY_WRITING(fp, name) \
	 {   if (OPEN_FOR_BINARY_WRITING(fp, name)) \
	     {   sprintf(error_msg, "opening '%s' for writing", name); \
		 return  CCPN_ERROR;   }   }

#define  CHECK_OPEN_FOR_WRITING(fp, name) \
	 {   if (OPEN_FOR_WRITING(fp, name)) \
	     {   sprintf(error_msg, "opening '%s' for writing", name); \
		 return  CCPN_ERROR;   }   }

#define  FCLOSE(fp)  {   if (fp)  fclose(fp);  fp = (FILE *) NULL;   }

#define  DISCARD_REST_OF_LINE(fp) \
	 {   Line L;  fgets(L, LINE_SIZE, fp);   }

#define  LSEEK_RELATIVE(fd, nwords) \
	 (((nwords) != 0) \
	   ? lseek(fd, (off_t) (nwords)*BYTES_PER_WORD, SEEK_CUR) : 0)

#define  LSEEK_ABSOLUTE(fd, nwords) \
	 lseek(fd, (off_t) (nwords)*BYTES_PER_WORD, SEEK_SET)

extern int fseek_relative_safe(FILE *fp, long int offset);

#define  FSEEK_RELATIVE_SIZE(fp, nwords, size) \
	 fseek(fp, (long int) ((nwords)*(size)), SEEK_CUR)

#ifdef WIN64
#define  FSEEK_RELATIVE_BYTES(fp, nbytes) \
         _fseeki64(fp, (__int64) (nbytes), SEEK_CUR)
#else /* not WIN32 */
#define  FSEEK_RELATIVE_BYTES(fp, nbytes) \
	 fseeko(fp, (off_t) (nbytes), SEEK_CUR)
#endif /* end if WIN32 */

#define  FSEEK_RELATIVE(fp, nwords) \
	 fseek(fp, (long int) ((nwords)*BYTES_PER_WORD), SEEK_CUR)
/*
	 fseek_relative_safe(fp, (long int) (nwords)*BYTES_PER_WORD))
*/
/* below no good because fseek does not like negative offsets
   more than 256*1024*1024 bytes in size */
/* switch back now */
/* below no good because of bug in fopen for update mode
	 (((nwords) != 0) \
	   ? fseek(fp, (long int) ((nwords)*BYTES_PER_WORD), SEEK_CUR) : 0)
*/

#define  FSEEK_ABSOLUTE(fp, nwords) \
	 fseek(fp, (long int) ((nwords)*BYTES_PER_WORD), SEEK_SET)

#ifdef WIN64
#define  FSEEK_ABSOLUTE_BYTES(fp, nbytes) \
        _fseeki64(fp, (__int64) (nbytes), SEEK_SET)
#else /* not WIN32 */
#define  FSEEK_ABSOLUTE_BYTES(fp, nbytes) \
	 fseek(fp, (off_t) (nbytes), SEEK_SET)	 
#endif /* end if WIN32 */

#ifdef  FREAD
#undef  FREAD
#endif

#ifdef  FWRITE
#undef  FWRITE
#endif

#define  FREAD(x, size, n, fp) \
	 (fread((void *) (x), (size), (n), (fp)) != (n))

#define  FWRITE(x, size, n, fp) \
	 (fwrite((void *) (x), (size), (n), (fp)) != (n))

#ifdef  MALLOC
#undef  MALLOC
#endif
#define  MALLOC(ptr, type, num) \
	 {   if ( ((ptr)=(type *) malloc((unsigned) ((num)*sizeof(type)))) \
	           == NULL )  return  CCPN_ERROR; \
             DEBUG_CODE_HEX_INT("malloc", ptr, num)   }

#define  MALLOC_NEW(ptr, type, num) \
	 {   if ( ((ptr)=(type *) malloc((unsigned) ((num)*sizeof(type)))) \
	           == NULL )  return  NULL; \
             DEBUG_CODE_HEX_INT("malloc_new", ptr, num)   }

#define  STRING_MALLOC_COPY(string1, string2) \
	 {   int Len = strlen(string2) + 1; \
	     MALLOC(string1, char, Len); \
	     strcpy(string1, string2);   }

#define  STRING_MALLOC_NEW(string1, string2) \
	 {   int Len = strlen(string2) + 1; \
	     MALLOC_NEW(string1, char, Len); \
	     strcpy(string1, string2);   }

#define  MALLOC2(ptr, type, num1, num2) \
	 {   int I; \
	     MALLOC(ptr, type *, num1); \
	     for (I = 0; I < num1; I++) \
		 MALLOC((ptr)[I], type, num2);   }

#define  MALLOC_ZERO(ptr, type, num) \
	 {   int I; \
	     if ( ((ptr)=(type *) malloc((unsigned) ((num)*sizeof(type)))) \
	           				== NULL )  return  CCPN_ERROR; \
	     for (I = 0; I < (num); I++)  (ptr)[I] = (type) NULL; \
             DEBUG_CODE_HEX_INT("malloc_zero", ptr, num)   }

#define  MALLOC2_ZERO(ptr, type, num1, num2) \
	 {   int J; \
	     MALLOC(ptr, type *, num1); \
	     for (J = 0; J < num1; J++) \
		 MALLOC_ZERO((ptr)[J], type, num2);   }

#ifdef  FREE
#undef  FREE
#endif
#define  FREE(ptr, type) \
	 {   DEBUG_CODE_HEX("free", ptr) \
             if ((ptr) != (type *) NULL) \
	     {   free((type *) (ptr));  (ptr) = (type *) NULL;   }   }

#define  FREE2(ptr, type, num) \
	 {   int I; \
	     if ((ptr) != (type **) NULL) \
	     {   for (I = 0; I < num; I++) \
	             FREE((ptr)[I], type); \
	     	 FREE(ptr, type *);   }   }

#ifdef  REALLOC
#undef  REALLOC
#endif
#define  REALLOC(ptr, type, num) \
	 {   type *Ptr; \
             DEBUG_CODE_HEX_INT("realloc1", ptr, num) \
	     if ( ((Ptr)=(type *) realloc((ptr), \
		    (unsigned)((num)*sizeof(type)))) \
	               == NULL )  return  CCPN_ERROR;  else  ptr = Ptr; \
             DEBUG_CODE_HEX("realloc2", ptr)   }

#define  REALLOC_NEW(ptr, type, num) \
	 {   type *Ptr; \
             DEBUG_CODE_HEX_INT("realloc_new1", ptr, num) \
	     if ( ((Ptr)=(type *) realloc((ptr), \
		    (unsigned)((num)*sizeof(type)))) \
	               == NULL )  return  NULL;  else  ptr = Ptr; \
             DEBUG_CODE_HEX("realloc_new2", ptr)   }

#define  CHECK_FOR_REALLOC(ptr, type, num, block) \
	 {   if ( !((num) % (block)) )  REALLOC(ptr, type, num+block);   }

#define  ZERO_VECTOR(v, n) \
	 {   int  I;  for (I = 0; I < (n); I++)  (v)[I] = 0;   }

#define  SUBTRACT_VECTORS(v1, v2, v3, n) \
	 {   int  I;  for (I = 0; I < (n); I++)  (v1)[I] = (v2)[I]-(v3)[I];   }

#define  ADD_VECTORS(v1, v2, v3, n) \
	 {   int  I;  for (I = 0; I < (n); I++)  (v1)[I] = (v2)[I]+(v3)[I];   }

#define  SCALE_VECTOR(v1, v2, s, n) \
	 {   int  I;  for (I = 0; I < (n); I++)  (v1)[I] = (s)*(v2)[I];   }

#define  INNER_PRODUCT(d, v1, v2, n) \
	 {   int  I;  d = 0; \
	     for (I = 0; I < (n); I++)  d += (v1)[I]*(v2)[I];   }

#define  VECTOR_PRODUCT(d, v, n) \
	 {   int  I;  d = 1; \
	     for (I = 0; I < (n); I++)  d *= (v)[I];   }

#define  CROSS_PRODUCT(v1, v2, v3) \
	 {   (v1)[0] = (v2)[1]*(v3)[2] - (v2)[2]*(v3)[1]; \
	     (v1)[1] = (v2)[2]*(v3)[0] - (v2)[0]*(v3)[2]; \
	     (v1)[2] = (v2)[0]*(v3)[1] - (v2)[1]*(v3)[0];   }

#define  NEGATE_VECTOR(v1, v2, n) \
	 {   int  I;  for (I = 0; I < (n); I++)  (v1)[I] = - (v2)[I];   }

#define  COPY_VECTOR(v1, v2, n) \
	 {   int  I;  for (I = 0; I < (n); I++)  (v1)[I] = (v2)[I];   }

#define  COPY_VECTOR_FROM_TOP(v1, v2, n) \
	 {   int  I;  for (I = (n)-1; I >= 0; I--)  (v1)[I] = (v2)[I];   }

#define  INDEX_OF_ARRAY(index, array, cumul, n) \
	 {   INNER_PRODUCT(index, array, cumul, n);  }

#define  ARRAY_OF_INDEX(array, index, cumul, n) \
	 {   int I; long Ind = index; \
	     for (I = (n)-1; I >= 0; I--) \
	     {   array[I] = Ind / cumul[I];  Ind %= cumul[I];   }   }

#define  CUMULATIVE(cumul, array, total, n) \
	 {   int I;  total = 1; \
	     for (I = 0; I < n; I++) \
	     {   (cumul)[I] = total;  total *= (array)[I];   }   }

#define  BLOCK(point, size)  ((int) (1 + ((point) - 1)/(size)))

#define  BLOCKS(blocks, points, size, n) \
	 {   int I; \
	     for (I = 0; I < n; I++) \
	     {   (blocks)[I] = BLOCK((points)[I], (size)[I]);   }   }

#define  STRIP_LEADING_SPACE(string) \
	 {   char *Ptr1 = string, *Ptr2 = string; \
	     while (*Ptr2 && isspace(*Ptr2))  Ptr2++; \
	     while (*Ptr2)  *Ptr1++ = *Ptr2++;  *Ptr1 = 0;   }

#define  STRIP_TRAILING_SPACE(string) \
	 {   if (strlen(string)) \
	     {   char *Ptr = string + strlen(string) - 1; \
		 while ((Ptr >= string) && isspace(*Ptr))  Ptr--; \
		 *++Ptr = 0;   };   }

#define  STRIP_SPACE(string) \
	{   STRIP_LEADING_SPACE(string);   STRIP_TRAILING_SPACE(string);   }

#define  STRIP_CARRIAGE_RETURN(string) \
	 {   char *Ptr; \
	     if ((Ptr = strchr(string, '\n')) != NULL)  *Ptr = 0;   }

/* use stdout instead of stderr so that Unix novices */
/* do not lose output when running jobs in background */
#define  ERROR_AND_EXIT(message) \
	 {   fprintf(stdout, "Fatal error: %s\n", message);  exit  (1);   }

#define  ERROR_AND_RETURN_ERROR(message) \
	 {   char Msg[200]; \
	     sprintf(Msg, "Error: %s", message); \
	     print_error_message(Msg);  return  CCPN_ERROR;   }

#define  ERROR_AND_RETURN(message) \
	 {   char Msg[200]; \
	     sprintf(Msg, "Error: %s", message); \
	     print_error_message(Msg);  return;   }

#define  RETURN_ERROR_MSG(message) \
	 {   sprintf(error_msg, message);  return  CCPN_ERROR;   }

#define  RETURN_ERROR_MSG_NULL(message) \
	 {   sprintf(error_msg, message);  return  NULL;   }

#define  RETURN_WITH_ERROR(msg, fmt, vars) \
	 {   sprintf(msg, fmt, vars);  return  CCPN_ERROR;   }

#define  ERROR_AND_ABORT(message) \
	 {   fprintf(stderr, "Error: %s\n", message);  abort();   }

#define  ASSERT(expression, message) \
	 {   if (!(expression)) \
	     {   char  full_message[LINE_SIZE]; \
		 sprintf(full_message, "Assertion '%s' failed\n", message); \
		 ERROR_AND_ABORT(full_message);   }   }

#define  CHECK_STATUS(status) \
	 {   if ((status) == CCPN_ERROR)  return  CCPN_ERROR;   }

#define  CHECK_NULL(null) \
	 {   if (!(null))  return NULL;   }

#define  CHECK_STATUS_NULL(status) \
	 {   if ((status) == CCPN_ERROR)  return  NULL;   }

#define  CHECK_NULL_STATUS(null) \
	 {   if (!(null))  return  CCPN_ERROR;   }

#define  CHECK_OK(status) \
	 {   CcpnStatus S;  if ((S = (status)) != CCPN_OK)  return  S;   }

#define  CHECK_RANGE(string, value, first, last) \
	 {   if (((value) < (first)) || ((value) > (last))) \
	     {   sprintf(error_msg, \
			"invalid choice for %s, must be in range %d-%d", \
						string, first, last); \
	     	 return  CCPN_ERROR;   }   }

#define  CHECK_LOWER(string, value, first) \
	 {   if ((value) < (first)) \
	     {   sprintf(error_msg, "invalid choice for %s, must be >= %d", \
							string, first); \
	     	 return  CCPN_ERROR;   }   }

#define  CHECK_UPPER(string, value, last) \
	 {   if ((value) > (last)) \
	     {   sprintf(error_msg, "invalid choice for %s, must be <= %d", \
							string, last); \
	     	 return  CCPN_ERROR;   }   }
/*
#define  CHECK_RANGE(string, value, first, last) \
	 {   char Msg[LINE_SIZE]; \
	     if (((value) < (first)) || ((value) > (last))) \
	     {   sprintf(Msg, "invalid choice for %s, must be in range %d-%d", \
						string, first, last); \
	     	 ERROR_AND_RETURN_ERROR(Msg);   }   }

#define  CHECK_LOWER(string, value, first) \
	 {   char Msg[LINE_SIZE]; \
	     if ((value) < (first)) \
	     {   sprintf(Msg, "invalid choice for %s, must be >= %d", \
							string, first); \
	     	 ERROR_AND_RETURN_ERROR(Msg);   }   }

#define  CHECK_UPPER(string, value, last) \
	 {   char Msg[LINE_SIZE]; \
	     if ((value) > (last)) \
	     {   sprintf(Msg, "invalid choice for %s, must be <= %d", \
							string, last); \
	     	 ERROR_AND_RETURN_ERROR(Msg);   }   }
*/

#endif /* _incl_macros */
