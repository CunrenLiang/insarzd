//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


//./simamp z.rdr sim.rdr 15070 3.0 100

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

FILE *openfile(char *filename, char *pattern);
void readdata(void *data, size_t blocksize, FILE *fp);
void writedata(void *data, size_t blocksize, FILE *fp);
float *array1d_float(long nc);
void free_array1d_float(float *fv);
double *array1d_double(long nc);
void free_array1d_double(double *fv);

int main(int argc, char *argv[]){

  FILE *hgtfp;
  FILE *simfp;

  double *hgt;
  float *sim;

  float scale;
  float offset;

  int nrg, naz;
  int i, j;


  if(argc != 6){
    fprintf(stderr, "\nUsage: %s hgt sim nrg scale offset\n\n", argv[0]);
    exit(1);
  }

  hgtfp = openfile(argv[1], "rb");
  simfp = openfile(argv[2], "wb");
  nrg = atoi(argv[3]);
  scale = atof(argv[4]);
  offset = atof(argv[5]);

  fseeko(hgtfp,0L,SEEK_END);
  naz = ftello(hgtfp) / sizeof(double) / nrg;
  rewind(hgtfp);

  printf("nrg: %d, naz: %d\n", nrg, naz);

  hgt = array1d_double(nrg);
  sim = array1d_float(nrg);

  for(i = 0; i < naz; i++){
    if((i + 1) % 100 == 0)
      fprintf(stderr,"processing line: %6d of %6d\r", i+1, naz);
    readdata((double *)hgt, nrg * sizeof(double), hgtfp);
    for(j = 0; j < nrg - 1; j++){
      sim[j] = (hgt[j+1] - hgt[j]) * scale + offset;
    }
    sim[nrg-1] = 0.0;
    writedata((float *)sim, nrg * sizeof(float), simfp);
  }
  fprintf(stderr,"processing line: %6d of %6d\n", naz, naz);

  free_array1d_double(hgt);
  free_array1d_float(sim);
  fclose(hgtfp);
  fclose(simfp);

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

double *array1d_double(long nc){

  double *fv;

  fv = (double*) malloc(nc * sizeof(double));
  if(!fv){
    fprintf(stderr,"Error: cannot allocate 1-D double vector\n");
    exit(1);
  }

  return fv;
}

void free_array1d_double(double *fv){
  free(fv);
}