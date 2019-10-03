//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


#include "resamp.h"

int main(int argc, char *argv[]){

  FILE *infp;   //slave image to be resampled
  FILE *outfp;  //resampled slave image

  fcomplex *filter;
  fcomplex *in;
  fcomplex *out;
  fcomplex *tmp;
  fcomplex *tmp2;
  int *zeroflag;

  int nrg; //file width
  int naz; //file length

  int nfft; //fft length
  int nfilter; //filter length
  int hnfilter;

  float bw;
  float bc;
  float beta; //kaiser window beta

  int zero_cf;
  float offset;

  float sc; //constant to scale the data read in to avoid large values
            //during fft and ifft
  float cf_pha;
  float t;
  fcomplex cf;

  int nblock_in;
  int nblock_out;
  int num_block;
  int i_block;
  int nblock_in_last;
  int nblock_out_last;

  int i, j;



/*****************************************************************************/
  nfilter = 65;
  nfft = 1024;
  beta = 1.0;
  zero_cf = 0;
  offset = 0.0;

  sc = 10000.0;
/*****************************************************************************/


  if(argc < 6){
    fprintf(stderr, "\nusage: %s inputf outputf nrg bw bc [nfilter] [beta] [zero_cf]\n", argv[0]);
    fprintf(stderr, "\nmandatory:\n");
    fprintf(stderr, "  inputf:  input file\n");
    fprintf(stderr, "  outputf: output file\n");
    fprintf(stderr, "  nrg:     file width\n");
    fprintf(stderr, "  bw:      filter bandwidth divided by sampling frequency [0, 1]\n");
    fprintf(stderr, "  bc:      filter center frequency divided by sampling frequency\n");
    fprintf(stderr, "\noptional:\n");
    fprintf(stderr, "  nfilter: number samples of the filter (odd). Default: 65\n");
    fprintf(stderr, "  nfft:    number of samples of the FFT. Default: 1024\n");
    fprintf(stderr, "  beta:    kaiser window beta. Default: 1.0\n");
    fprintf(stderr, "  zero_cf: if bc != 0.0, move center frequency to zero? 0: Yes (Default). 1: No.\n");
    fprintf(stderr, "  offset:  offset (in samples) of linear phase for moving center frequency. Default: 0.0\n\n");
    exit(1);
  }

  //open files
  infp  = openfile(argv[1], "rb");
  outfp = openfile(argv[2], "wb");

  nrg  = atoi(argv[3]);
  naz  = file_length(infp, nrg, sizeof(fcomplex));
  printf("file width: %d, file length: %d\n\n", nrg, naz);

  bw = atof(argv[4]);
  bc = atof(argv[5]);
  
  if(argc > 6)
    nfilter = atoi(argv[6]);
  if(argc > 7)
    nfft = atoi(argv[7]);

  if(argc > 8)
    beta = atof(argv[8]);
  if(argc > 9)
    zero_cf = atoi(argv[9]);
  if(argc > 10)
    offset = atof(argv[10]);

  if(nfilter < 3){
    fprintf(stderr, "filter length: %d too small!\n", nfilter);
    exit(1);
  }
  if(nfilter % 2 != 1){
    fprintf(stderr, "filter length must be odd!\n");
    exit(1);
  }

  hnfilter = (nfilter - 1) / 2;

  nblock_in = nfft - nfilter + 1;

  nblock_in += hnfilter;
  if (nblock_in <= 0){
    fprintf(stderr, "fft length too small compared with filter length!\n");
    exit(1);
  }

  nblock_out = nblock_in - 2 * hnfilter;

  num_block = (nrg - 2 * hnfilter) / nblock_out;
  if((nrg - num_block * nblock_out - 2 * hnfilter) != 0){
    num_block += 1;
  }
  if((nrg - 2 * hnfilter) <= 0){
    num_block = 1;
  }
  if(num_block == 1){
    nblock_out_last = 0;
    nblock_in_last = nrg;
  }
  else{
    nblock_out_last = nrg - (num_block - 1) * nblock_out - 2 * hnfilter;
    nblock_in_last = nblock_out_last + 2 * hnfilter;
  }

  filter = array1d_fcomplex(nfft);
  in     = array1d_fcomplex(nrg);
  out    = array1d_fcomplex(nrg);
  tmp    = array1d_fcomplex(nfft);
  tmp2   = array1d_fcomplex(nfft);
  zeroflag = array1d_int(nrg);


  bandpass_filter(bw, bc, nfilter, nfft, (nfilter-1)/2, beta, filter);

  //relationship of nr and matlab fft
  //nr fft           matlab fft
  //  1      <==>     ifft()*nfft
  // -1      <==>     fft()


  four1((float *)filter - 1, nfft, -1);


  for(i = 0; i < naz; i++){

    if((i + 1) % 1000 == 0 || (i + 1) == naz)
      fprintf(stderr,"processing line: %6d of %6d\r", i+1, naz);
    if((i + 1) == naz)
      fprintf(stderr,"\n\n");
  
    //read data
    readdata((fcomplex *)in, (size_t)nrg * sizeof(fcomplex), infp);
    for(j = 0; j < nrg; j++){
      if(in[j].re != 0.0 || in[j].im != 0.0){
        zeroflag[j] = 1;
      }
      else{
        zeroflag[j] = 0;
      }
      in[j].re *= 1.0 / sc;
      in[j].im *= 1.0 / sc;
    }

    //process
    for(i_block = 0; i_block < num_block; i_block++){
      //zero out
      for(j = 0; j < nfft; j++){
        tmp[j].re = 0.0;
        tmp[j].im = 0.0;
      }

      //get data
      if(num_block == 1){
        for(j = 0; j < nrg; j++){
          tmp[j] = in[j];
        }
      }
      else{
        if(i_block == num_block - 1){
          for(j = 0; j < nblock_in_last; j++){
            tmp[j] = in[j+nblock_out*i_block];
          }
        }
        else{
          for(j = 0; j < nblock_in; j++){
            tmp[j] = in[j+nblock_out*i_block];
          }
        }
      }


      four1((float *)tmp - 1, nfft, -1);

      for(j = 0; j < nfft; j++)
        tmp2[j] = cmul(filter[j], tmp[j]);

      four1((float *)tmp2 - 1, nfft, 1);


      if(num_block == 1){
        for(j = 0; j < nrg; j++){
          out[j] = tmp2[j];
        }
      }
      else{
        if(i_block == 0){
          for(j = 0; j < hnfilter + nblock_out; j++){
            out[j] = tmp2[j];
          }
        }
        else if(i_block == num_block - 1){
          for(j = 0; j < hnfilter + nblock_out_last; j++){
            out[nrg - 1 - j] = tmp2[nblock_in_last - 1 - j];
          }
        }
        else{
          for(j = 0; j < nblock_out; j++){
            out[j + hnfilter + i_block * nblock_out] = tmp2[j + hnfilter];
          }
        }
      }
    }

    //write data
    //move center frequency
    if(bc != 0 && zero_cf == 0){
      for(j = 0; j < nrg; j++){
        //t = j - (nrg - 1.0) / 2.0; //make 0 index exactly at range center
        t = j + offset; //make 0 index exactly at range center
        cf_pha = 2.0 * PI * (-bc) * t;
        cf.re = cos(cf_pha);
        cf.im = sin(cf_pha);
        out[j] = cmul(out[j], cf);
      }
    }
    //scale back
    for(j = 0; j < nrg; j++){
      if(zeroflag[j] == 0){
        out[j].re *= 0.0;
        out[j].im *= 0.0;
      }
      else{
        out[j].re *= sc / nfft;
        out[j].im *= sc / nfft;     
      }
    }
    //write data
    writedata((fcomplex *)out, nrg * sizeof(fcomplex), outfp);
  }

  free_array1d_fcomplex(filter);
  free_array1d_fcomplex(in);
  free_array1d_fcomplex(out);
  free_array1d_fcomplex(tmp);
  free_array1d_fcomplex(tmp2);
  free_array1d_int(zeroflag);
  fclose(infp);
  fclose(outfp);

  return 0;
}//end main()


