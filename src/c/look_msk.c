//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


#include "resamp.h"

int main(int argc, char *argv[]){

  FILE *infp;
  FILE *outfp;
  signed char **in;
  double *a;
  double sum;
  signed char *out;
  long nrg, naz;
  long nrg1, naz1;
  int nrlks, nalks;
  int i, j, k;

  if(argc < 4){
    fprintf(stderr, "\nUsage: %s infile outfile nrg nrlks nalks\n\n", argv[0]);
    fprintf(stderr, "  infile:  input file\n");
    fprintf(stderr, "  outfile: output file\n");
    fprintf(stderr, "  nrg:     file width\n");
    fprintf(stderr, "  nrlks:   number of looks in range (default: 4)\n");
    fprintf(stderr, "  nalks:   number of looks in azimuth (default: 4)\n\n");
    exit(1);
  }

  infp  = openfile(argv[1], "rb");
  outfp = openfile(argv[2], "wb");
  
  nrg = atoi(argv[3]);
  naz = file_length(infp, nrg, sizeof(signed char));

  if(argc > 4)
    nrlks = atoi(argv[4]);
  else
    nrlks = 4;

  if(argc > 5)
    nalks = atoi(argv[5]);
  else
    nalks = 4;
  
  nrg1 = nrg / nrlks;
  naz1 = naz / nalks;

  in = array2d_char(nalks, nrg);
  a  = array1d_double(nrg);
  out= array1d_char(nrg1);


  for(i = 0; i < naz1; i++){

    if((i + 1) % 100 == 0)
      fprintf(stderr,"processing line: %6d of %6d\r", i+1, naz1);

    readdata((signed char *)in[0], (size_t)nalks * (size_t)nrg * sizeof(signed char), infp);
    //take looks in azimuth
    for(j = 0; j < nrg; j++){
      a[j] = 0.0;
      for(k = 0; k < nalks; k++){
        a[j] += in[k][j];
      }
    }
    //take looks in range
    for(j = 0; j < nrg1; j++){
      sum = 0.0;
      for(k = 0; k < nrlks; k++){
        sum += a[j * nrlks + k];
      }
      out[j] = (signed char)(sum / nrlks / nalks);
    }
  
    writedata((signed char *)out, nrg1 * sizeof(signed char), outfp);
  }
  fprintf(stderr,"processing line: %6d of %6d\n", naz1, naz1);

  free_array2d_char(in);
  free_array1d_double(a);
  free_array1d_char(out);
  fclose(infp);
  fclose(outfp);

  return 0;
}
