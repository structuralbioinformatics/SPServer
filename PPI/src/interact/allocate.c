#include     "interact.h"

 void nrerror(error_text)
 char error_text[];
 {
 void exit();
 fprintf(stderr,"Numerical Recipes run-time error...\n");
 fprintf(stderr,"%s\n",error_text);
 fprintf(stderr,"...now exiting to system...\n");
 exit(1);
 }

float  **Fmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 float  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(float **)  malloc((unsigned) (nrh-nrl+1)*sizeof(float *));
 if (!m) nrerror("allocation failure 1 in Fmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
  if (!m[i]) nrerror("allocation failure 2 in Fmatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_Fmatrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((float *) (m[i]+ncl));
   free((float *) (m+nrl));
}

