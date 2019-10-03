//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


// program for correcting phase
// Programmer: Cunren Liang

// This program calculates the average differential phase of two consecutive frames,
// then fits a polynomial to the differential phase, and finally corrects the second
// frame's phase with the fitted differential phase.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NORMAL_EXIT 0
#define ABNORMAL_EXIT 0

#define NR_END 1
#define FREE_ARG char*
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

typedef struct {
  float re;
  float im;
} fcomplex;

typedef struct {
  double re;
  double im;
} dcomplex;

//
FILE *openfile(char *filename, char *pattern);
void readdata(void *data, size_t blocksize, FILE *fp);
void writedata(void *data, size_t blocksize, FILE *fp);
long file_length(FILE* fp, long width, long element_size);
fcomplex *array1d_fcomplex(long nc);
void free_array1d_fcomplex(fcomplex *fcv);
dcomplex *array1d_dcomplex(long nc);
void free_array1d_dcomplex(dcomplex *fcv);
float *array1d_float(long nc);
void free_array1d_float(float *fv);
float *vector_float(long nl, long nh);
void free_vector_float(float *v, long nl, long nh);
fcomplex cmul(fcomplex a, fcomplex b);
fcomplex cconj(fcomplex z);
//
void covsrt(float **covar, int ma, int ia[], int mfit);
void gaussj(float **a, int n, float **b, int m);
void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
  int ma, float **covar, float *chisq, void (*funcs)(float, float [], int));
void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
//
void funcs(float x,float afunc[],int ma);
int *vector_int(long nl, long nh);
void free_vector_int(int *v, long nl, long nh);
float **matrix_float(long nrl, long nrh, long ncl, long nch);
void free_matrix_float(float **m, long nrl, long nrh, long ncl, long nch);


int main(int argc, char *argv[]){

  FILE *infp0, *infp1;
  FILE *outfp;

  fcomplex *in0, *in1;
  fcomplex *out;

  int nrgin0, nrgin1;
  int nazin0, nazin1;
  int nrgout;
  int nazout;

  int nrgoff;
  int nazoff;

  int rgs0, rge0; //starting and ending range indexes (included) of overlap area in image 0
  int rgs1, rge1; //starting and ending range indexes (included) of overlap area in image 1
  int azs0, aze0; //starting and ending azimuth indexes (included) of overlap area in image 0
  int azs1, aze1;  //starting and ending azimuth indexes (included) of overlap area in image 1
  int nrgo; //range length of overlap area in image 0
  int nazo; //azimuth length of overlap area in image 1

  dcomplex *diff; //differential phase of overlap area

  int ds, de; //starting and ending indexes (included) of diff for polynomial fitting
  int nd; // number of data points for fitting

  //variables for polynomial fitting
  float *x;  // containning portion of C style indices [0...nrgin1-1]
  float *y;  // differential phase
  float *sig; //stardard deviations of data points in x and y

  float *a; //polynomial coefficients
  int *ia;  //indicates whether this polynomial coefficents should be fitted or not
  int ma; // number of polynomial coefficients

  float **covar; //covariance matrix of fitting
  float chisq;

  float *afunc; //values of polynomial base functions
  float *phafit; //phase calculated from fitted polynomial

  int index;
  fcomplex tmp;
  double mag;
  double mag1, mag2;

  fcomplex ax, bx;
  dcomplex cx;

  int i, j;

  int order; //order of polynomial
  order = 3;
  ma = order + 1;


  if(argc != 8){
    fprintf(stderr, "\nUsage: %s input0 nrgin0 input1 nrgin1 nrgoff nazoff output\n\n", argv[0]);
    fprintf(stderr, "input0: first image (e.g. upper frame)\n");
    fprintf(stderr, "nrgin0: width of inputfile0\n");
    fprintf(stderr, "input1: second image (e.g. lower frame)\n");
    fprintf(stderr, "nrgin1: width of inputfile1\n");
    fprintf(stderr, "nrgoff: range offset in samples (input0 is master)\n");
    fprintf(stderr, "nazoff: azimuth offset in samples (input0 is master)\n");
    fprintf(stderr, "output: processed second image with phase corrected by phase difference of the two images\n\n");
    exit(1);
  }

  //read inputs:
  infp0  = openfile(argv[1], "rb");
  nrgin0 = atoi(argv[2]);
  nazin0 = file_length(infp0, nrgin0, sizeof(fcomplex));

  infp1  = openfile(argv[3], "rb");
  nrgin1 = atoi(argv[4]);
  nazin1 = file_length(infp1, nrgin1, sizeof(fcomplex));

  nrgoff = atoi(argv[5]);
  nazoff = atoi(argv[6]);
  outfp = openfile(argv[7], "wb");

  //check parameters:
  if(nrgoff > nrgin1 - 1 || nrgin0 - 1 + nrgoff < 0){
    fprintf(stderr, "Error: There is no overlap in range direction\n\n");
    exit(1);
  }

  if(nazoff > nazin1 - 1 || nazin0 - 1 + nazoff < 0){
    fprintf(stderr, "Error: There is no overlap in azimuth direction\n\n");
    exit(1);
  }

  in0  = array1d_fcomplex(nrgin0);
  in1  = array1d_fcomplex(nrgin1);
  out  = array1d_fcomplex(nrgin1);
  diff = array1d_dcomplex(nrgin1);
  phafit = array1d_float(nrgin1);

  ////////////////////////////////////////////////////////////////////////
  // STEP 1. find starting and ending indices for differential processing
  ////////////////////////////////////////////////////////////////////////

  //finding range starting index
  for(i = 0; i <= nrgin1 - 1; i++){
    index = i - nrgoff;
    if(index >= 0 && index <= nrgin0 - 1){
      rgs0 = index;
      rgs1 = i;
      break;
    }
  }
  //finding range ending index
  for(i = nrgin1 - 1; i >=0; i--){
    index = i - nrgoff;
    if(index >= 0 && index <= nrgin0 - 1){
      rge0 = index;
      rge1 = i;
      break;
    }
  }
  //finding azimuth starting index
  for(i = 0; i <= nazin1 - 1; i++){
    index = i - nazoff;
    if(index >= 0 && index <= nazin0 - 1){
      azs0 = index;
      azs1 = i;
      break;
    }
  }
  //finding azimuth ending index
  for(i = nazin1 - 1; i >= 0; i--){
    index = i - nazoff;
    if(index >= 0 && index <= nazin0 - 1){
      aze0 = index;
      aze1 = i;
      break;
    }
  }

  nrgo = rge1 - rgs1 + 1;
  nazo = aze1 - azs1 + 1;



  ////////////////////////////////////////////////////////////////////////
  // STEP 2. calculate average differential phase
  ////////////////////////////////////////////////////////////////////////

  fprintf(stderr,"calculate average phase difference\n");

  for(i = 0; i < nrgin1; i++){
    diff[i].re = 0.0;
    diff[i].im = 0.0;
  }

  //skip non-overlap area
  fseeko(infp0, (size_t)azs0 * (size_t)nrgin0 * sizeof(fcomplex), SEEK_SET);
  fseeko(infp1, (size_t)azs1 * (size_t)nrgin1 * sizeof(fcomplex), SEEK_SET);

  //calculate phase difference
  for(i = 0; i < nazo; i++){

    if((i + 1) % 100 == 0)
      fprintf(stderr,"processing line: %6d of %6d\r", i + 1, nazo);
    if(i + 1 == nazo)
      fprintf(stderr,"processing line: %6d of %6d\n\n", i + 1, nazo);

    readdata((fcomplex *)in1, nrgin1 * sizeof(fcomplex), infp1);
    readdata((fcomplex *)in0, nrgin0 * sizeof(fcomplex), infp0);
    for(j = rgs1; j < rge1; j++){
      if(in0[j - nrgoff].re != 0.0 && in0[j - nrgoff].im != 0.0 && in1[j].re != 0.0 && in1[j].im != 0.0){
        
        //nomalize to avoid too large values which may cause errors
        //CL, AUG-2015
        //not a good idea. comment out. CL, SEP-2015.
        //mag1 = sqrt(in0[j - nrgoff].re * in0[j - nrgoff].re + in0[j - nrgoff].im * in0[j - nrgoff].im);
        //mag2 = sqrt(in1[j].re * in1[j].re + in1[j].im * in1[j].im);
        //in1[j].re /= mag2;
        //in1[j].im /= mag2;
        //in0[j - nrgoff].re /= mag1;
        //in0[j - nrgoff].im /= mag1;

        //tmp = cmul(in0[j - nrgoff], cconj(in1[j]));

        //remove one of the magintude to improve interferogram averaging. CL JUN-2016
        mag2 = sqrt(in1[j].re * in1[j].re + in1[j].im * in1[j].im);
        
        in1[j].im = -in1[j].im;

        ax = in0[j - nrgoff];
        bx = in1[j];
        cx.re=ax.re*bx.re-ax.im*bx.im;
        cx.im=ax.im*bx.re+ax.re*bx.im;

        diff[j].re += cx.re / mag2;
        diff[j].im += cx.im / mag2;
      }
    }
  }

  //normalize
  for(i = 0; i < nrgin1; i++){
    if(diff[i].re != 0.0 || diff[i].im != 0.0){
      mag = sqrt(diff[i].re * diff[i].re + diff[i].im * diff[i].im);
      diff[i].re /= mag;
      diff[i].im /= mag;
    }
  }

  //find phase for fitting
  for(i = 0; i <= nrgin1 - 1; i++){
    if(diff[i].re != 0.0 && diff[i].im != 0.0){
      ds = i;
      break;
    }
  }
  for(i = nrgin1 - 1; i >= 0; i--){
    if(diff[i].re != 0.0 && diff[i].im != 0.0){
      de = i;
      break;
    }
  }
  nd = de - ds + 1;



  ////////////////////////////////////////////////////////////////////////
  // STEP 3. fit a polynomial to the average differential phase
  ////////////////////////////////////////////////////////////////////////

  fprintf(stderr,"fit a polynomial to the differential phase\n\n");

  //prepare phase for fitting
  x = vector_float(1, nd);
  y = vector_float(1, nd);
  for(i = ds; i <= de; i++){
    x[i - ds + 1] = i;
    y[i - ds + 1] = atan2(diff[i].im, diff[i].re);
  }

  //setting equal standard deviations
  sig = vector_float(1, nd);
  for(i = 1; i <= nd; i++){
    sig[i] = 1.0;
  }

  //setting polynomial coefficents to be fitted for
  a   = vector_float(1, ma);
  ia  = vector_int(1, ma);
  for(i = 1; i <= ma; i++){
    ia[i] = 1;
  }
  
  //allocate covariance matrix
  covar = matrix_float(1, ma, 1, ma);
  
  //fit
  lfit(x, y, sig, nd, a, ia, ma, covar, &chisq, funcs);



  ////////////////////////////////////////////////////////////////////////
  // STEP 4. correct phase of next image using fitted phase
  ////////////////////////////////////////////////////////////////////////

  fprintf(stderr,"correct phase of next frame\n");

  //calculate fitted differential phase
  afunc = vector_float(1, ma);
  for(i = 0; i < nrgin1; i++){
    //if(i >= ds && i <= de){
      phafit[i] = 0.0;
      funcs(x[i - ds + 1], afunc, ma);
      for(j = 1; j <= ma; j++){
        phafit[i] += a[j] * afunc[j];
      }
    //}
    //else{
    //  phafit[i] = 0.0;
    //}
  }

  //correct next image's phase using fitted phase
  rewind(infp1);
  for(i = 0; i < nazin1; i++){

    if((i + 1) % 1000 == 0)
      fprintf(stderr,"processing line: %6d of %6d\r", i + 1, nazin1);
    if(i + 1 == nazin1)
      fprintf(stderr,"processing line: %6d of %6d\n\n", i + 1, nazin1);

    readdata((fcomplex *)in1, nrgin1 * sizeof(fcomplex), infp1);
    for(j = 0; j < nrgin1; j++){
      //if(phafit[j] != 0.0 && in1[j].re != 0.0 && in1[j].im != 0.0){
      if(in1[j].re != 0.0 && in1[j].im != 0.0){
        tmp.re = cos(phafit[j]);
        tmp.im = sin(phafit[j]);
        out[j] = cmul(in1[j], tmp);
      }
      else{
        out[j].re = 0.0;
        out[j].im = 0.0;
      }
    }
    writedata((fcomplex *)out, nrgin1 * sizeof(fcomplex), outfp);
  }


  free_array1d_fcomplex(in0);
  free_array1d_fcomplex(in1);
  free_array1d_fcomplex(out);

  free_array1d_dcomplex(diff);
  free_array1d_float(phafit);

  free_vector_float(x, 1, nd);
  free_vector_float(y, 1, nd);
  free_vector_float(sig, 1, nd);

  free_vector_float(a, 1, ma);
  free_vector_int(ia, 1, ma);
  free_matrix_float(covar, 1, ma, 1, ma);
  free_vector_float(afunc, 1, ma);

  fclose(infp0);
  fclose(infp1);
  fclose(outfp);

  return 0;

}

////////////////////////////////////////////////////////////////////////////
FILE *openfile(char *filename, char *pattern){
  FILE *fp;
  
  fp=fopen(filename, pattern);
  if (fp==NULL){
    fprintf(stderr,"Error: cannot open file: %s\n", filename);
    exit(1);
  }

  return fp;
}

void readdata(void *data, size_t blocksize, FILE *fp){
  if(fread(data, blocksize, 1, fp) != 1){
    fprintf(stderr,"Error: cannot read data\n");
    exit(1);
  }
}

void writedata(void *data, size_t blocksize, FILE *fp){
  if(fwrite(data, blocksize, 1, fp) != 1){
    fprintf(stderr,"Error: cannot write data\n");
    exit(1);
  }
}

long file_length(FILE* fp, long width, long element_size){
  long length;
  
  fseeko(fp,0L,SEEK_END);
  length = ftello(fp) / element_size / width;
  rewind(fp);
  
  return length;
}

fcomplex *array1d_fcomplex(long nc){

  fcomplex *fcv;

  fcv = (fcomplex*) malloc(nc * sizeof(fcomplex));
  if(!fcv){
    fprintf(stderr,"Error: cannot allocate 1-D float vector\n");
    exit(1);
  }

  return fcv;
}

void free_array1d_fcomplex(fcomplex *fcv){
  free(fcv);
}

dcomplex *array1d_dcomplex(long nc){

  dcomplex *fcv;

  fcv = (dcomplex*) malloc(nc * sizeof(dcomplex));
  if(!fcv){
    fprintf(stderr,"Error: cannot allocate 1-D double complex vector\n");
    exit(1);
  }

  return fcv;

}

void free_array1d_dcomplex(dcomplex *fcv){
  free(fcv);
}

float *array1d_float(long nc){

  float *fv;

  fv = (float*) malloc(nc * sizeof(float));
  if(!fv){
    fprintf(stderr,"Error: cannot allocate 1-D float vector\n");
    exit(1);
  }

  return fv;
}

void free_array1d_float(float *fv){
  free(fv);
}

float *vector_float(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v){
    fprintf(stderr,"Error: cannot allocate 1-D vector\n");
    exit(1);  
  }
  
  return v-nl+NR_END;
}

void free_vector_float(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

fcomplex cmul(fcomplex a, fcomplex b)
{
  fcomplex c;
  c.re=a.re*b.re-a.im*b.im;
  c.im=a.im*b.re+a.re*b.im;
  return c;
}

fcomplex cconj(fcomplex z)
{
  fcomplex c;
  c.re=z.re;
  c.im = -z.im;
  return c;
}
////////////////////////////////////////////////////////////////////////////

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void covsrt(float **covar, int ma, int ia[], int mfit)
{
  int i,j,k;
  float swap;

  for (i=mfit+1;i<=ma;i++)
    for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit;
  for (j=ma;j>=1;j--) {
    if (ia[j]) {
      for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
      for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
      k--;
    }
  }
}
#undef SWAP

#define NRANSI
//#include "nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(float **a, int n, float **b, int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  float big,dum,pivinv,temp;

  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for (j=1;j<=n;j++) ipiv[j]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j] != 1)
        for (k=1;k<=n;k++) {
          if (ipiv[k] == 0) {
            if (fabs(a[j][k]) >= big) {
              big=fabs(a[j][k]);
              irow=j;
              icol=k;
            }
          }
        }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
      for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
        dum=a[ll][icol];
        a[ll][icol]=0.0;
        for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
        for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
        SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI

#define NRANSI
//#include "nrutil.h"

void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
  int ma, float **covar, float *chisq, void (*funcs)(float, float [], int))
{
  void covsrt(float **covar, int ma, int ia[], int mfit);
  void gaussj(float **a, int n, float **b, int m);
  int i,j,k,l,m,mfit=0;
  float ym,wt,sum,sig2i,**beta,*afunc;

  beta=matrix(1,ma,1,1);
  afunc=vector(1,ma);
  for (j=1;j<=ma;j++)
    if (ia[j]) mfit++;
  if (mfit == 0) nrerror("lfit: no parameters to be fitted");
  for (j=1;j<=mfit;j++) {
    for (k=1;k<=mfit;k++) covar[j][k]=0.0;
    beta[j][1]=0.0;
  }
  for (i=1;i<=ndat;i++) {
    (*funcs)(x[i],afunc,ma);
    ym=y[i];
    if (mfit < ma) {
      for (j=1;j<=ma;j++)
        if (!ia[j]) ym -= a[j]*afunc[j];
    }
    sig2i=1.0/SQR(sig[i]);
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
        wt=afunc[l]*sig2i;
        for (j++,k=0,m=1;m<=l;m++)
          if (ia[m]) covar[j][++k] += wt*afunc[m];
        beta[j][1] += ym*wt;
      }
    }
  }
  for (j=2;j<=mfit;j++)
    for (k=1;k<j;k++)
      covar[k][j]=covar[j][k];
  gaussj(covar,mfit,beta,1);
  for (j=0,l=1;l<=ma;l++)
    if (ia[l]) a[l]=beta[++j][1];
  *chisq=0.0;
  for (i=1;i<=ndat;i++) {
    (*funcs)(x[i],afunc,ma);
    for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
    *chisq += SQR((y[i]-sum)/sig[i]);
  }
  covsrt(covar,ma,ia,mfit);
  free_vector(afunc,1,ma);
  free_matrix(beta,1,ma,1,1);
}
#undef NRANSI

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) nrerror("cannot allocate vector()");
  return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("cannot allocate ivector()");
  return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("cannot allocate matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) nrerror("cannot allocate matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


////////////////////////////////////////////////////////////////////////////
void funcs(float x,float afunc[],int ma)
{
  //This is a polynomial
  // a1 + a2 * x + a3 * x^2 + a4 * x^3 +... 
  int i;

  for(i = 1; i <= ma; i++){
    afunc[i] = pow(x, i - 1.0);
  }
}

int *vector_int(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("Error: cannot allocate vector_int()");
  return v-nl+NR_END;
}

void free_vector_int(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

float **matrix_float(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("Error: cannot allocate vector2d_float()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) nrerror("Error: cannot allocate vector2d_float()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_matrix_float(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

  

