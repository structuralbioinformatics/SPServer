#include <stdio.h>
#include <math.h>

#define TRUE        ((unsigned char) 0xff)

#define	NR_END	1
#define	SIGN(a, b)	((b) >= 0.0 ? fabs(a) : -fabs(a))

static	double	sqrarg;
#define SQR(a)	((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)


void	tred2(),  tqli();
void	nrerror();
void	shell();
int	*ivector();
double	**matrix(), *vector();
double	pythag(double a, double b);

main (argc, argv)
int	argc;
char	*argv[];
{
	int	i, j; 

	FILE	*outfile;

	a = matrix(1, (long)natm, 1, (long)natm);
	d = vector(1, (long)natm);
	d_ent = vector(1, (long)natm);
	e = vector(1, (long)natm);
	ind = ivector(1, (long)natm); 
	atlist = ivector(1, (long)natm);

           tred2(a, nwnatm, d, e);
           tqli(d, e, nwnatm, a);

/* Sorting the eigenvalues
*/

	shell(nwnatm, d, ind);

	fprintf(outfile, "   VAR.     EIG.VAL.        FREQ(s-1)   FREQ.(cm-1)  %c VARIAB.     CUM.VARIAB.   CONT.ENTR.     CUMUL.ENTR.       T.DS (Kcal/mol)      REM. T.DS (Kcal/mol)\n\n", percnt);
	var_par = 0.;
	cum_par = 0.;
        cum_ent = 0.;
	for (i = nwnatm; i >= nmin; i--) {
	    cum_ent += d_ent[ind[i]];
            var_par = 100. * d[ind[i]] / tot_eig;
	    cum_par += var_par;
	    if (d[ind[i]] <= LowFreq) {
	       freqs = -1.;
               freqcm = -1.;
	    }
	    else {
               freqs = sqrt(KT / (d[ind[i]] * AM));
               freqcm = freqs / CLIGHT;
            }
	    fprintf(outfile, " %5d   %14.7f    %10.3e   %10.3f     %6.2f         %6.2f %14.7f %14.7f    %14.7f             %14.7f\n", nwnatm - i + 1, d[ind[i]], freqs, freqcm, var_par, cum_par, d_ent[ind[i]], cum_ent, cum_ent * Temp, (tot_ent - cum_ent) * Temp);
	}
	fclose(outfile);

/* Writing the eigenvectors
*/

	   ivar = 0;
	   outfile = fopen(ei_vc_fnm, "w");
	   fprintf(outfile, " ");
	   for (i = 1; i <= nwnatm; i++) {
	       for (j = 1; j <= nwnatm; j++) {
                   fprintf(outfile, " %14.7E", a[j][ind[i]]); 
		   ivar++;
	           if (ivar == 5) {
	   	      ivar = 0;
	              fprintf(outfile, "\n ");
	           }
	       }
	   }  
	   fclose(outfile);

}
/***********************************/
void	tred2(a, n, d, e)
double	**a, *d, *e;
int	n;

{
	int	l, k, j, i;
	double	scale, hh, h, g, f;

	for (i = n; i >= 2; i--) {
	    l = i - 1;
	   
	    h = scale = 0.0;

	    if (l > 1) {

		  for (k = 1; k <= l; k++) {
			 scale += fabs(a[i][k]);
		  }

		  if (scale == 0.0) {
			e[i] = a[i][l];
		  }
		  else {
			for (k = 1; k <= l; k++) {
			    a[i][k] /= scale;
			    h += a[i][k] * a[i][k];
			}
			f = a[i][l];
			g = (f >= 0.0 ? -sqrt(h) : sqrt(h));

			e[i] = scale * g;

			h -= f * g;

			a[i][l] = f - g;

			f = 0.0;

			for (j = 1; j <= l; j++) {
			    a[j][i] = a[i][j] / h;
			    g = 0.0;

			    for (k = 1; k <= j; k++) {
				   g += a[j][k] * a[i][k];
			    }

			    for (k = j + 1; k <= l; k++) {
				   g += a[k][j] * a[i][k];
			    }
	
			    e[j] = g / h;

			    f += e[j] * a[i][j];

			}

			hh = f / (h + h);

			for (j = 1; j <= l; j++) {
			    f = a[i][j];
			    e[j] = g = e[j] - hh * f;

                            for (k = 1; k <= j; k++) {
    	                        a[j][k] -= (f * e[k] + g * a[i][k]);
			    }
			}

		  }

	    }

	    else {
		    e[i] = a[i][l];
	    }

	    d[i] = h;

	}

	d[1] = 0.0;
	e[1] = 0.0;

	for (i = 1; i <= n; i++) {
	    l = i - 1;
	    if (d[i]) {

	       for (j = 1; j <= l; j++) {
	   	   g = 0.0;
		   for (k = 1; k <= l; k++) {
		       g += a[i][k] * a[k][j];
		   }

		   for (k = 1; k <= l; k++) {
		       a[k][j] -= g * a[k][i];
		   }

	       }
	    }

	    d[i] = a[i][i];
		 
            a[i][i] = 1.0;

	    for (j = 1; j <= l; j++) {
	        a[j][i] = a[i][j] = 0.0;
	    }

	}

}
/***********************************/
void tqli(double *d, double *e, int n, double **z)
{
	int	m, l, iter, i, k;	
	double	s, r, p, g, f, dd, c, b;

	for (i = 2; i <= n; i++) {
	    e[i - 1] = e[i];
	}

	e[n] = 0.0;

	for (l = 1; l <= n; l++) {
	    iter = 0;

	    do {
	       for (m = l; m <= n - 1; m++) {
		   dd = fabs(d[m]) + fabs(d[m + 1]);
		   if ((double)(fabs(e[m]) + dd) == dd) break;
	       }

	       if (m != l) {
		  if (iter++ == 30) nrerror("Too many iterations in tqli");
		  g = (d[l + 1] - d[l]) / (2.0 * e[l]);

		  r = pythag(g, 1.0);
		  g = d[m] - d[l] + e[l] / (g + SIGN(r, g));

		  s = c = 1.0;
		  p = 0.0; 

		  for (i = m - 1; i >= l; i--) {
		      f = s * e[i];
		      b = c * e[i];

		      e[i + 1] = (r = pythag(f, g));

		      if (r == 0.0) {
			 d[i + 1] -= p;
			 e[m] = 0.0;
			 break;
		      }

		      s = f / r;
		      c = g / r;

		      g = d[i + 1] - p;
		      r = (d[i] - g) * s + 2.0 * c * b;

		      d[i + 1] = g + (p =s*r);
		      g = c * r - b;

		      for (k = 1; k <= n; k++) {
			  f = z[k][i + 1];
			  z[k][i + 1] = s * z[k][i] + c * f;
			  z[k][i] = c * z[k][i] - s * f;
		      }
	          }

	          if (r == 0.0 && i >= l) continue;
	          d[l] -= p;
	          e[l] = g;
	          e[m] = 0.0;

	       }

	    } while (m != l);
	}
}
/***********************************/
double **matrix(long nrl, long nrh, long ncl, long nch) 
{
	long	i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double	**m;

	m = (double **) malloc((size_t)((nrow + NR_END)*sizeof(double*)));

	if (!m) nrerror("Allocation failure 1 in matrix()");

	m += NR_END;
	m -= nrl;

	m[nrl] = (double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));


	if (!m[nrl]) nrerror("Allocation failure 2 in matrix()"); 
	
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) {
	    m[i] = m[i - 1] + ncol;
	}

	return m;
	
}
/***********************************/
double *vector(long nl, long nh)
{

	double	*v;

	v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END)*sizeof(double)));

	if (!v) nrerror("Allocation failure in vector()");

	return v - nl + NR_END;
}
/***********************************/
int	*ivector(long nl, long nh)
{

        int	*v;

        v = (int *) malloc((size_t) ((nh - nl + 1 + NR_END)*sizeof(int)));
        if (!v) nrerror("Allocation failure in vector()");

        return v - nl + NR_END;
}
/***********************************/
void nrerror(char error_text[])
{
	
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}
/***********************************/
double pythag(double a, double b)
{
	double	absa, absb;

	absa = fabs(a);
	absb = fabs(b);

	if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}
/*-- Sorting subroutine called shell, taken from Numerical Recipes in C, Press et al., pp 244-245, and modified by X.de la Cruz --*/
void shell(n, arr, ind)
int     n, ind[];
double	arr[];
{
        int     i, j, m, nn, lognb2, it;
        double	t;

        lognb2 = (log((double) n) * ALN2I + TINY);
        m = n;
        for (i = 1; i <= n; i++) ind[i] = i;

        for (nn = 1; nn <= lognb2; nn++) {
            m >>= 1;
            for (j = m + 1; j <= n; j++) {
                i = j - m;
                it = ind[j];
                t = arr[it];
                while (i >= 1 && arr[ind[i]] > t) {
                      ind[i + m] = ind[i];
                      i -= m;
                }
                ind[i + m] = it;
            }
        }
}
/***********************************/
 char  *cvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   char  *v;
   void nrerror();
   v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
   if (!v) nrerror("allocation failure in cvector");
   return v-nl;
 }

void free_cvector(v,nl,nh)
char  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((char*) (v+nl));
}


