//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct {
  float re;
  float im;
} fcomplex;

FILE *openfile(char *filename, char *pattern);
void readdata(void *data, size_t blocksize, FILE *fp);
void writedata(void *data, size_t blocksize, FILE *fp);
fcomplex cmul(fcomplex a, fcomplex b);
fcomplex cconj(fcomplex z);
float xcabs(fcomplex z);
fcomplex *array1d_fcomplex(long nc);
void free_array1d_fcomplex(fcomplex *fcv);

int main(int argc, char *argv[]){

  FILE *slc1fp;
  FILE *slc2fp;
  FILE *intfp;
  FILE *ampfp;

  fcomplex *slc1;
  fcomplex *slc2;
  fcomplex *intf;
  fcomplex *amp;

  int nrg;
  int naz1;
  int naz2;

  int i,j;

  if(argc != 6){
    fprintf(stderr, "\nUsage: %s slc1 slc2 intf amp nrg\n\n", argv[0]);
    exit(1);
  }

  slc1fp = openfile(argv[1], "rb");
  slc2fp = openfile(argv[2], "rb");
  intfp = openfile(argv[3], "wb");
  ampfp = openfile(argv[4], "wb");
  nrg = atoi(argv[5]);

  fseeko(slc1fp,0L,SEEK_END);
  naz1 = ftello(slc1fp) / sizeof(fcomplex) / nrg;
  fseeko(slc2fp,0L,SEEK_END);
  naz2 = ftello(slc2fp) / sizeof(fcomplex) / nrg;
  rewind(slc1fp);
  rewind(slc2fp);

  if(naz1 != naz2){
    fprintf(stderr, "WARNING: file sizes are different\n");
  }
  printf("width: %d, length: %d\n", nrg, naz1);

  slc1 = array1d_fcomplex(nrg);
  slc2 = array1d_fcomplex(nrg);
  intf = array1d_fcomplex(nrg);
  amp  = array1d_fcomplex(nrg);

  for(i = 0; i < naz1; i++){

    if((i + 1) % 100 == 0)
      fprintf(stderr,"processing line: %6d of %6d\r", i+1, naz1);

    readdata((fcomplex *)slc1, nrg * sizeof(fcomplex), slc1fp);
    readdata((fcomplex *)slc2, nrg * sizeof(fcomplex), slc2fp);

    for(j = 0; j < nrg; j++){
      intf[j] = cmul(slc1[j], cconj(slc2[j]));
      amp[j].re = xcabs(slc1[j]);
      amp[j].im = xcabs(slc2[j]);
    }

    writedata((fcomplex *)intf, nrg * sizeof(fcomplex), intfp);
    writedata((fcomplex*)amp, nrg * sizeof(fcomplex), ampfp);
  }
  fprintf(stderr,"processing line: %6d of %6d\n", naz1, naz1);

  free_array1d_fcomplex(slc1);
  free_array1d_fcomplex(slc2);
  free_array1d_fcomplex(intf);
  free_array1d_fcomplex(amp);
  fclose(slc1fp);
  fclose(slc2fp);
  fclose(intfp);
  fclose(ampfp);

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

float xcabs(fcomplex z)
{
  float x,y,ans,temp;
  x=fabs(z.re);
  y=fabs(z.im);
  if (x == 0.0)
    ans=y;
  else if (y == 0.0)
    ans=x;
  else if (x > y) {
    temp=y/x;
    ans=x*sqrt(1.0+temp*temp);
  } else {
    temp=x/y;
    ans=y*sqrt(1.0+temp*temp);
  }
  return ans;
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







