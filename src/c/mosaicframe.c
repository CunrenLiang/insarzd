//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


// program for mosaicing multiple consecutive frames
// Cunren Liang, 22-MAY-2015
// JPL/Caltech


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

FILE **array1d_FILE(long nc);
void free_array1d_FILE(FILE **fv);
int *array1d_int(long nc);
void free_array1d_int(int *fv);
fcomplex *array1d_fcomplex(long nc);
void free_array1d_fcomplex(fcomplex *fcv);
fcomplex **array2d_fcomplex(long nl, long nc);
void free_array2d_fcomplex(fcomplex **m);

long file_length(FILE* fp, long width, long element_size);

int main(int argc, char *argv[]){

  FILE **infp;
  FILE *outfp;
  
  fcomplex *in;
  fcomplex *out, *out1, *out2;  
  fcomplex **upperoverlap; //upper overlap area;
  fcomplex **loweroverlap; //lower overlap area;

  int n;
  int *nrgin; //infile width
  int *nazin; //infile length
  int *nrgoff; //infile range offset
  int *nazoff; //infile azimuth offset
  int *oflag; //overlap output flag
  int nrgout; //outfile width
  int nazout; //outfile length

  int uos, uoe, uol; //start, end and length of upper overlap area
  int los, loe, lol; //start, end and length of lower overlap area
  int cnl; //length of center area

  int paroff;
  int parcyc;

  char diffname[256];
  FILE *difffp;
  fcomplex *diff;
  int diffflag;
  diffflag = 0;

  int i, j, k;

  int delta; //edge to be removed of the overlap area (number of lines)
  delta = 20;

  if(argc < 5){
    fprintf(stderr, "\nUsage: %s outputfile nrgout nazout n [inputfile0] [nrgin0] [nrgoff0] [nazoff0] [oflag0] (repeat...)\n\n", argv[0]);
    fprintf(stderr, "  for each frame\n");
    fprintf(stderr, "  range offset is relative to the left most sample\n");
    fprintf(stderr, "  azimuth offset is relative to the first line of last frame\n\n");
    exit(1);
  }

  //read mandatory parameters
  outfp  = openfile(argv[1], "wb");
  nrgout = atoi(argv[2]);
  nazout = atoi(argv[3]);
  n      = atoi(argv[4]);

  //allocate memory
  infp   = array1d_FILE(n);
  nrgin  = array1d_int(n);
  nazin  = array1d_int(n);
  nrgoff = array1d_int(n); //nrgoff must be <= 0
  nazoff = array1d_int(n); //nazoff must be <= 0
  oflag  = array1d_int(n);

  in   = array1d_fcomplex(nrgout);
  out  = array1d_fcomplex(nrgout);
  out1 = array1d_fcomplex(nrgout);
  out2 = array1d_fcomplex(nrgout);
  diff = array1d_fcomplex(nrgout);

  //read optional parameters
  paroff = 4;
  parcyc = 5;
  for(i = 0; i < n; i++){
    infp[i] = openfile(argv[paroff + parcyc*i + 1], "rb");
    nrgin[i] = atoi(argv[paroff + parcyc*i + 2]);
    nrgoff[i] = atoi(argv[paroff + parcyc*i + 3]);
    nazoff[i] = atoi(argv[paroff + parcyc*i + 4]);
    oflag[i] = atoi(argv[paroff + parcyc*i + 5]);
    nazin[i] = file_length(infp[i], nrgin[i], sizeof(fcomplex));
    if(nrgoff[i] > 0){
      fprintf(stderr,"Error: positive range offset: %d\n\n", nrgoff[i]);
      exit(1);
    }
    if(nazoff[i] > 0){
      fprintf(stderr,"Error: positive azimuth offset: %d\n\n", nazoff[i]);
      exit(1);
    }
    if(nrgout < nrgin[i] - nrgoff[i]){
      fprintf(stderr,"Error: ouput length < nrgin[%d] - nrgoff[%d], %d, %d\n\n", i, i, nrgout, nrgin[i] - nrgoff[i]);
      exit(1);
    }
  }


  for(i = 0; i < n; i++){

    fprintf(stderr,"processing frame: %2d of %2d\n", i+1, n);
    
    //////////////////////////////////////////////////////////////////////
    // STEP 1. calculate overlap area
    //         split the whole area into three parts: upper overlap area, 
    //         center non-overlap area, and lower overlap area
    //////////////////////////////////////////////////////////////////////

    //we follow the following convention: line and column number start with 0.
    //upper overlap area of frame i
    if(i != 0){
      uos = - nazoff[i];
      uoe = nazin[i-1] - 1;
      uol = uoe - uos + 1;
      if(uol < delta * 2){
        fprintf(stderr,"Error: not enough overlap area between frame: %d and %d\n\n", i-1, i);
        exit(1);
      }
    }
    else{
      uos = 0;
      uoe = 0;
      uol = 0;
    }
    
    //lower overlap area of frame i
    if(i != n - 1){
      los = - nazoff[i+1];
      loe = nazin[i] - 1;
      lol = loe - los + 1;
      if(lol < delta * 2){
        fprintf(stderr,"Error: not enough overlap area between frame: %d and %d\n\n", i, i+1);
        exit(1);
      }
    }
    else{
      los = 0;
      loe = 0;
      lol = 0;
    }
    
    //center non-overlap area of frame i
    //should add a check here?
    cnl = nazin[i] - uol - lol;


    //////////////////////////////////////////////////////////////////////
    // STEP 2. write out center non-overlap area of frame i
    //         for frame 0 (first frame),
    //         this is the area only exluding the lower overlap area
    //         for last frame
    //         this is the area only excluding the upper overlap area
    //////////////////////////////////////////////////////////////////////

    //prepare for reading and writing data
    for(k = 0; k < nrgout; k++){
      in[k].re = 0.0;
      in[k].im = 0.0;
      out[k].re = 0.0;
      out[k].im = 0.0;
    }

    //read and write center non-overlap area
    //this only excludes the lower overlap area for the first frame
    //this only excludes the upper overlap area for the last frame
    for(j = 0; j < cnl; j++){
      readdata((fcomplex *)in, nrgin[i] * sizeof(fcomplex), infp[i]);
      for(k = 0; k < nrgin[i]; k++){
        out[k - nrgoff[i]].re = in[k].re;
        out[k - nrgoff[i]].im = in[k].im;
      }
      writedata((fcomplex *)out, nrgout * sizeof(fcomplex), outfp);
    }


    //////////////////////////////////////////////////////////////////////
    // STEP 3. deal with lower overlap area of frame i
    //         which corresponds to upper overlap area of frame i + 1
    //////////////////////////////////////////////////////////////////////

    //deal with last frame
    if(i == n - 1){
    //  for(j = 0; j < lol; j++){
    //    for(k = 0; k < nrgin[i]; k++){
    //      out[k - nrgoff[i]].re = loweroverlap[j][k].re;
    //      out[k - nrgoff[i]].im = loweroverlap[j][k].im;
    //    } 
    //    writedata((fcomplex *)out, nrgout * sizeof(fcomplex), outfp);
    //  }
    //  free_array2d_fcomplex(loweroverlap);
      
      //no need to do anything
      break;
    }

    sprintf(diffname, "%d-%d.int", i, i+1);
    difffp = openfile(diffname, "wb");

    loweroverlap = array2d_fcomplex(lol, nrgin[i]);
    readdata((fcomplex *)loweroverlap[0], lol * nrgin[i] * sizeof(fcomplex), infp[i]);
    upperoverlap = array2d_fcomplex(lol, nrgin[i+1]);
    readdata((fcomplex *)upperoverlap[0], lol * nrgin[i+1] * sizeof(fcomplex), infp[i+1]);

    //deal with other frames
    for(j = 0; j < lol; j++){

      //prepare for writing data
      for(k = 0; k < nrgout; k++){
        out[k].re = 0.0;
        out[k].im = 0.0;
        out1[k].re = 0.0;
        out1[k].im = 0.0;
        out2[k].re = 0.0;
        out2[k].im = 0.0;
        diff[k].re = 0.0;
        diff[k].im = 0.0;
      }

      //upper edge of overlap area
      //use data of current frame: frame i
      if(j < delta){
        for(k = 0; k < nrgin[i]; k++){
          out[k - nrgoff[i]].re = loweroverlap[j][k].re;
          out[k - nrgoff[i]].im = loweroverlap[j][k].im;
        }
      }
      //lower edge of overlap area
      //use data of next frame: frame i+1
      else if(j >= lol - delta){
        for(k = 0; k < nrgin[i+1]; k++){
          out[k - nrgoff[i+1]].re = upperoverlap[j][k].re;
          out[k - nrgoff[i+1]].im = upperoverlap[j][k].im;
        }
      }
      //center area of overlap area
      else{
        
        //prepare data
        for(k = 0; k < nrgin[i]; k++){
          out1[k - nrgoff[i]].re = loweroverlap[j][k].re;
          out1[k - nrgoff[i]].im = loweroverlap[j][k].im;
        }
        for(k = 0; k < nrgin[i+1]; k++){
          out2[k - nrgoff[i+1]].re = upperoverlap[j][k].re;
          out2[k - nrgoff[i+1]].im = upperoverlap[j][k].im;
        }

        //output difference of overlap area
        //diffflag 0: upper frame phase - lower frame phase
        if(diffflag == 0){
          for(k = 0; k < nrgout; k++)
            if(out1[k].re != 0 && out1[k].im != 0 && out2[k].re != 0 && out2[k].im != 0){
              diff[k] = cmul(out1[k], cconj(out2[k]));
            }
        }
        //diffflag 1: upper frame - lower frame
        else {
          for(k = 0; k < nrgout; k++)
            if(out1[k].re != 0 && out1[k].im != 0 && out2[k].re != 0 && out2[k].im != 0){
              diff[k].re = out1[k].re - out2[k].re;
              diff[k].im = out1[k].im - out2[k].im;
            }
        }
        writedata((fcomplex *)diff, nrgout * sizeof(fcomplex), difffp);
        
        //mosaic overlap area
        //case0: mosaic at the center of overlap area
        if(oflag[i] == 0){
          if(j < lol/2){
            for(k = 0; k < nrgout; k++){
              out[k].re = out1[k].re;
              out[k].im = out1[k].im;
            }
          }
          else{
            for(k = 0; k < nrgout; k++){
              out[k].re = out2[k].re;
              out[k].im = out2[k].im;
            }
          }
        }
        //case1: mosaic at the lower edge of overlap area
        else if(oflag[i] == 1){
          for(k = 0; k < nrgout; k++){
            out[k].re = out1[k].re;
            out[k].im = out1[k].im;
          }
        }
        //case2: mosaic at the upper edge of overlap area
        else if(oflag[i] == 2){
          for(k = 0; k < nrgout; k++){
            out[k].re = out2[k].re;
            out[k].im = out2[k].im;
          }
        }
        //case3: average
        else{
          for(k = 0; k < nrgout; k++){
            out[k].re = out1[k].re + out2[k].re;
            out[k].im = out1[k].im + out2[k].im;
            if(out1[k].re && out1[k].im && out2[k].re && out2[k].im){
              out[k].re /= 2.0;
              out[k].im /= 2.0;
            }
          }
        }
      }
      writedata((fcomplex *)out, nrgout * sizeof(fcomplex), outfp);
    }

    free_array2d_fcomplex(loweroverlap);
    free_array2d_fcomplex(upperoverlap);
  
    fclose(difffp);
  }

  for(i = 0; i < n; i++)
    fclose(infp[i]);
  fclose(outfp);

  free_array1d_FILE(infp);
  free_array1d_int(nrgin);
  free_array1d_int(nazin);
  free_array1d_int(nrgoff); //nrgoff must be <= 0
  free_array1d_int(nazoff); //nazoff must be <= 0
  free_array1d_int(oflag);

  free_array1d_fcomplex(in);
  free_array1d_fcomplex(out);
  free_array1d_fcomplex(out1);
  free_array1d_fcomplex(out2);

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

FILE **array1d_FILE(long nc){

  FILE **fv;

  fv = (FILE **)malloc(nc * sizeof(FILE *));
  if(!fv){
    fprintf(stderr,"Error: cannot allocate 1-D FILE array\n");
    exit(1);
  }

  return fv;
}

void free_array1d_FILE(FILE **fv){
  free(fv);
}

int *array1d_int(long nc){

  int *fv;

  fv = (int*) malloc(nc * sizeof(int));
  if(!fv){
    fprintf(stderr,"Error: cannot allocate 1-D int array\n");
    exit(1);
  }

  return fv;
}

void free_array1d_int(int *fv){
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

fcomplex **array2d_fcomplex(long nl, long nc){
/* allocate a fcomplex 2-D matrix */

  fcomplex **m;
  int i;

  /* allocate pointers to rows */
  m = (fcomplex **) malloc(nl * sizeof(fcomplex *));
  if(!m){
    fprintf(stderr,"Error: cannot allocate 2-D matrix\n");
    exit(1);
  }
 
  /* allocate rows */ 
  m[0] = (fcomplex*) malloc(nl * nc * sizeof(fcomplex));
  if(!m[0]){
    fprintf(stderr,"Error: cannot allocate 2-D matrix\n");
    exit(1);
  }

   /* set pointers */
  for(i = 1; i < nl; i++){
    m[i] = m[i-1] + nc;
  }

  return m;
}

void free_array2d_fcomplex(fcomplex **m){
/* free a fcomplex matrix allocated by fcarray2d() */
  free(m[0]);
  free(m);
}

long file_length(FILE* fp, long width, long element_size){
  long length;
  
  fseeko(fp,0L,SEEK_END);
  length = ftello(fp) / element_size / width;
  rewind(fp);
  
  return length;
}
