//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI 3.1415926535897932384626433832795028841971693993751058

typedef struct {
  float re;
  float im;
} fcomplex;


FILE *openfile(char *filename, char *pattern);
void readdata(void *data, size_t blocksize, FILE *fp);
void writedata(void *data, size_t blocksize, FILE *fp);
float *array1d_float(long nc);
void free_array1d_float(float *fv);
fcomplex *array1d_fcomplex(long nc);
void free_array1d_fcomplex(fcomplex *fcv);
fcomplex cmul(fcomplex a, fcomplex b);
fcomplex cconj(fcomplex z);

int main(int argc, char *argv[]){

  FILE *intfp;
  FILE *rgoffp;
  FILE *flatfp;

  fcomplex *intf;
  float *rgoff;
  fcomplex *flat;

  float rgs;
  float wvl;
  float flat_phs;
  fcomplex tmp;

  long nsf;
  long nsf2;
  long nsb;
  long nsr;
  long nb;

  int i, j;

  
  nsb = 1024 * 1024 * 10; //10M

  if(argc != 6){
    fprintf(stderr, "\nUsage: %s inteferogram rangeoffset flatterned_interferogram rgs wvl\n\n", argv[0]);
    exit(1);
  }

  intfp = openfile(argv[1], "rb");
  rgoffp = openfile(argv[2], "rb");
  flatfp = openfile(argv[3], "wb");
  fseeko(intfp,0L,SEEK_END);
  nsf = ftello(intfp) / sizeof(fcomplex);
  rewind(intfp);

  //check file size
  fseeko(rgoffp,0L,SEEK_END);
  nsf2 = ftello(rgoffp) / sizeof(float);
  rewind(rgoffp);
  if(nsf != nsf2){
    fprintf(stderr, "Error: file dimensions are different");
    exit(1);
  }

  rgs = atof(argv[4]);
  wvl = atof(argv[5]);

  nb = nsf / nsb;
  nsr = nsf - nb * nsb;

  intf = array1d_fcomplex(nsb);
  rgoff = array1d_float(nsb);
  flat = array1d_fcomplex(nsb);
  

  for(i = 0; i < nb; i++){
    readdata((fcomplex *)intf, nsb * sizeof(fcomplex), intfp);
    readdata((float *)rgoff, nsb * sizeof(float), rgoffp);
    for(j = 0; j < nsb; j++){
      if(isnan(rgoff[j]) != 0 || 
         isnan(intf[j].re) != 0 || 
         isnan(intf[j].im) != 0 ||
         rgoff[j] == 0. || (intf[j].re == 0. && intf[j].im == 0.)){
        flat[j].re = 0.0;
        flat[j].im = 0.0;
      }
      else{
        flat_phs = -rgoff[j] * 4.0 * PI * rgs / wvl;
        tmp.re = cos(flat_phs);
        tmp.im = sin(flat_phs);
        flat[j] = cmul(intf[j], tmp);
      }
    }
    writedata((fcomplex *)flat, nsb * sizeof(fcomplex), flatfp);
  }

  //dealing with last block
  if(nsr != 0){
    readdata((fcomplex *)intf, nsr * sizeof(fcomplex), intfp);
    readdata((float *)rgoff, nsr * sizeof(float), rgoffp);
    for(j = 0; j < nsr; j++){
      if(isnan(rgoff[j]) != 0 || 
         isnan(intf[j].re) != 0 || 
         isnan(intf[j].im) != 0 ||
         rgoff[j] == 0. || (intf[j].re == 0. && intf[j].im == 0.)){
        flat[j].re = 0.0;
        flat[j].im = 0.0;
      }
      else{
        flat_phs = -rgoff[j] * 4.0 * PI * rgs / wvl;
        tmp.re = cos(flat_phs);
        tmp.im = sin(flat_phs);
        flat[j] = cmul(intf[j], tmp);
      }
    }
    writedata((fcomplex *)flat, nsr * sizeof(fcomplex), flatfp);
  }

  free_array1d_fcomplex(intf);
  free_array1d_fcomplex(flat);
  free_array1d_float(rgoff);
  fclose(intfp);
  fclose(flatfp);
  fclose(rgoffp);

  return 0;
}

/******************************************************************/
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

double *array1d_double(long nc){

  double *fv;

  fv = (double*) malloc(nc * sizeof(double));
  if(!fv){
    fprintf(stderr,"Error: cannot allocate 1-D float vector\n");
    exit(1);
  }

  return fv;
}

void free_array1d_double(double *fv){
  free(fv);
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