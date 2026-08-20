/*
======================COPYRIGHT/LICENSE START==========================

diag_dbl.c: Part of the CcpNmr Analysis program

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
#include "diag_dbl.h"

static void tred2(double **a, int n,double *d, double *e)
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	d--; e--;
        for (i = 0; i < n; i++)
		a[i]--;
        a--;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = f>0 ? -sqrt(h) : sqrt(h);
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
				/* Next statement can be omitted if eigenvectors not wanted */
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	/* Next statement can be omitted if eigenvectors not wanted */
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}

	a++;
        for (i = 0; i < n; i++)
		a[i]++;
}

#define FSIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
    
static int tqli(double *d, double *e, int n, double **z)
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;
	int status = 1;
    
	d--; e--;
	for (i = 0; i < n; i++)
		z[i]--;
	z--;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if (m != l) {
				if (iter++ == 30)
				{
					status = 0;
					goto done;
				}

				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=sqrt((g*g)+1.0);
				g=d[m]-d[l]+e[l]/(g+FSIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					if (fabs(f) >= fabs(g)) {
						c=g/f;
						r=sqrt((c*c)+1.0);
						e[i+1]=f*r;
						c *= (s=1.0/r);
					} else {
						s=f/g;
						r=sqrt((s*s)+1.0);
						e[i+1]=g*r;
						s *= (c=1.0/r);
					}
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					p=s*r;
					d[i+1]=g+p;
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted */
					for (k=1;k<=n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				d[l]=d[l]-p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}

done:
	z++;
	for (i = 0; i < n; i++)
		z[i]++;

	return status;
}

CcpnStatus diagonalise_dbl(int n, double **matrix,
			double *eigenvalues, double *work, CcpnString error_msg)
{
    tred2(matrix, n, eigenvalues, work);

    if (!tqli(eigenvalues, work, n, matrix))
	RETURN_ERROR_MSG("diagonalisation did not converge");

    return  CCPN_OK;
}
