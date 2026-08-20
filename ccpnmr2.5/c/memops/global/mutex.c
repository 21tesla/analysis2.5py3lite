
/*
======================COPYRIGHT/LICENSE START==========================

mutex.c: Part of the CcpNmr Analysis program

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
/* implementation using posix threads */

/* this needs to go first because on IRIX this includes sys/param.h
   and that defines MIN and MAX without checking whether already defined */
#ifndef WIN64
#include <pthread.h>
#endif

#include "mutex.h"

#ifdef WIN64

#include <windows.h>

Mutex new_mutex()
{
    CRITICAL_SECTION *mutex;

    MALLOC_NEW(mutex, CRITICAL_SECTION, 1);

    InitializeCriticalSection(mutex);

    return (Mutex) mutex;
}

void delete_mutex(Mutex mutex)
{
    CRITICAL_SECTION *m = (CRITICAL_SECTION *) mutex;

    if (m)
    {
	DeleteCriticalSection(m);
	FREE(m, CRITICAL_SECTION);
    }
}

CcpnStatus lock_mutex(Mutex mutex)
{
    CRITICAL_SECTION *m = (CRITICAL_SECTION *) mutex;

    EnterCriticalSection(m);

    return CCPN_OK;
}

CcpnStatus unlock_mutex(Mutex mutex)
{
    CRITICAL_SECTION *m = (CRITICAL_SECTION *) mutex;

    LeaveCriticalSection(m);

    return CCPN_OK;
}

#else /* WIN32 */

Mutex new_mutex()
{
    pthread_mutex_t *mutex;

    MALLOC_NEW(mutex, pthread_mutex_t, 1);

    if (pthread_mutex_init(mutex, NULL))
    {
	FREE(mutex, pthread_mutex_t);

	return  (Mutex) NULL;
    }

    return (Mutex) mutex;
}

void delete_mutex(Mutex mutex)
{
    pthread_mutex_t *m = (pthread_mutex_t *) mutex;

    if (m)
    {
	(void) pthread_mutex_destroy(m);
	FREE(m, pthread_mutex_t);
    }
}

CcpnStatus lock_mutex(Mutex mutex)
{
    pthread_mutex_t *m = (pthread_mutex_t *) mutex;

    if (pthread_mutex_lock(m))
	return CCPN_ERROR;

    return CCPN_OK;
}

CcpnStatus unlock_mutex(Mutex mutex)
{
    pthread_mutex_t *m = (pthread_mutex_t *) mutex;

    if (pthread_mutex_unlock(m))
	return CCPN_ERROR;

    return CCPN_OK;
}

#endif /* WIN32 */

