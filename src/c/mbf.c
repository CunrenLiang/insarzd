//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


#include "resamp.h"

int main(int argc, char *argv[]){

  FILE *infp;
  FILE *outfp;

  fcomplex **in; //data read in
  fcomplex **out; //data written to output file
  fcomplex *filter; //multi-band bandpass filter
  fcomplex *filter_j;
  fcomplex *deramp; //deramp signal
  fcomplex *reramp; //reramp signal
  fcomplex *data; //data to be filtered.

  int nrg; //file width
  int naz; //file length
  float prf; //assume prf are the same
  float prf_frac; // azimuth processed bandwidth = prf_frac * prf
  float nb; //burst length in terms of pri. number of lines
  float nbg; //burst gap length in terms of pri. number of lines
  float nbc; //burst cycle length in terms of pri. number of lines
  float nboff; //number of unsynchronized lines in a burst with sign
              //see burst_sync.py for rules of sign.
              //the image to be processed is always considered to be master
              //and the other image is always considered to be slave
  float bsl; //burst start line, input float
  float kacoeff[3]; //FM rate along range (experessed by quadratic polynomial
                    //as a function of range sample number)
  float dopcoeff1[4]; //doppler centroid frequency along range (expressed by quadratic polynomial
                      //as a function of range sample number). this image
  float dopcoeff2[4];  //doppler centroid frequency along range (expressed by quadratic polynomial
                      //as a function of range sample number). the other image
                      //ATTENTION: MAKE RANGE NUMBER THE SAME ACCORDING RANGE OFFSET!!!

  float pri; // 1.0/prf
  float *ka;
  float *dop1;
  float *dop2;
  
  float *nfa; //full aperture length in terms of pri. number of lines
  float *freqs; //burst starting doppler centroid frequency
  float *freqe; //burst ending doppler centroid frequency
  float *bis; //burst imaged area start line numbers
  float *bie; //burst imaged area ending line numbers
  float *bic; //burst imaged area center line number, corresponding to the center of raw burst,
              //rather than the actual center of imaged area
  float *bica; //burst imaged area center line number, corresponding to the actual center of imaged area

  float deramp_center; //line number where center frequency is zero Hz after deramping

  float bis_min;
  float bis_max;
  float bie_min;
  float bie_max;

  int bis_out; //starting line number of the data block written out
  int bie_out; //ending line number of the data block written out
  int bis_in; //start line number of the data block read in
  int bie_in; //ending line number of the data block read in

  int bis_out2; //starting line number of the data block written out
  int bie_out2; //ending line number of the data block written out
  int bis_in2; //start line number of the data block read in
  int bie_in2; //ending line number of the data block read in

  float nb_new;
  float nbg_new;
  float nbc_new;
  float bsl_new;
  int nbc_new_int;

  int nburst_new; //number of bursts in a burst cycle

  float bfw; //bandwidth of burst in Hz
  float bfc; //center frequency of burst in Hz

  int nfft; //fft length
  int nfilter; //filter length, MUST BE ODD
  int hnfilter; //half filter length
  int edgl; //number of lines on the starting and ending edges

  float beta; //kaiser window beta
  float sc; //constant to scale the data read in to avoid large values
            //during fft and ifft
  int edgl_flag; //flag to indicate how many lines to keep on the starting and ending edges
                 //0: do not remove data on the edges
                 //1: remove data less than half convolution
                 //2: remove all data of incomplete convolution
  int deramp_center_flag; //flag to indicate the location with zero center frequency after
                          //deramping
                          //0: center (raw burst center) of the burst whose ending/start line number is used
                          //1: center of the burst cycle being processed
                          //2: center (raw burst center) of the center burst in the burst cycle being processed

  float tmp1, tmp2, tmp3;
  int i, j, k;

/*****************************************************************************/
//I just put these parametes which can be set here. These can also be set via
//arguments before running the programs if modifying the code to accept these
//arguments.

  beta = 1.0;
  nfilter = 257; //MUST BE ODD
  sc = 10000.0;
  edgl_flag = 0;
  deramp_center_flag = 0;
/*****************************************************************************/


  if(argc != 21){
    fprintf(stderr, "\nusage: %s inputf outputf nrg prf prf_frac nb nbg nboff bsl kacoeff[0-2] dopcoeff1[0-3] dopcoeff2[0-3]\n", argv[0]);
    fprintf(stderr, "\nmandatory:\n");
    fprintf(stderr, "  inputf:          input file\n");
    fprintf(stderr, "  outputf:         output file\n");
    fprintf(stderr, "  nrg:             file width\n");
    fprintf(stderr, "  prf:             PRF\n");
    fprintf(stderr, "  prf_frac:        fraction of PRF processed\n");
    fprintf(stderr, "                      (represents azimuth bandwidth)\n");
    fprintf(stderr, "  nb:              number of lines in a burst\n");
    fprintf(stderr, "                      (float, in terms of 1/PRF)\n");
    fprintf(stderr, "  nbg:             number of lines in a burst gap\n");
    fprintf(stderr, "                      (float, in terms of 1/PRF)\n");
    fprintf(stderr, "  nboff:           number of unsynchronized lines in a burst\n");
    fprintf(stderr, "                      (float, in terms of 1/PRF, with sign, see burst_sync.py for rules of sign)\n");
    fprintf(stderr, "                      (the image to be processed is always considered to be master)\n");
    fprintf(stderr, "  bsl:             start line number of a burst\n");
    fprintf(stderr, "                      (float, the line number of the first line of the full-aperture SLC is zero)\n");
    fprintf(stderr, "                      (no need to be first burst, any one is OK)\n");

    fprintf(stderr, "  kacoeff[0-2]:    FM rate coefficients\n");
    fprintf(stderr, "                      (three coefficients of a quadratic polynomial with regard to)\n");
    fprintf(stderr, "                      (range sample number. range sample number starts with zero)\n");

    fprintf(stderr, "  dopcoeff1[0-3]:  Doppler centroid frequency coefficients of this image\n");
    fprintf(stderr, "                      (four coefficients of a third order polynomial with regard to)\n");
    fprintf(stderr, "                      (range sample number. range sample number starts with zero)\n");

    fprintf(stderr, "  dopcoeff2[0-3]:  Doppler centroid frequency coefficients of the other image\n");
    fprintf(stderr, "                      (four coefficients of a third order polynomial with regard to)\n");
    fprintf(stderr, "                      (range sample number. range sample number starts with zero)\n\n");

    exit(1);
  }


  //open files
  infp  = openfile(argv[1], "rb");
  outfp = openfile(argv[2], "wb");

  nrg  = atoi(argv[3]);

  prf = atof(argv[4]);
  prf_frac = atof(argv[5]);
  nb = atof(argv[6]);
  nbg = atof(argv[7]);
  nboff = atof(argv[8]);
  bsl = atof(argv[9]);

  kacoeff[0] = atof(argv[10]);
  kacoeff[1] = atof(argv[11]);
  kacoeff[2] = atof(argv[12]);

  dopcoeff1[0] = atof(argv[13]);
  dopcoeff1[1] = atof(argv[14]);
  dopcoeff1[2] = atof(argv[15]);
  dopcoeff1[3] = atof(argv[16]);

  dopcoeff2[0] = atof(argv[17]);
  dopcoeff2[1] = atof(argv[18]);
  dopcoeff2[2] = atof(argv[19]);
  dopcoeff2[3] = atof(argv[20]);

  printf("\n\ninput parameters:\n");
  printf("input file: %s\n", argv[1]);
  printf("output file: %s\n", argv[2]);
  printf("nrg: %d\n", nrg);
  printf("prf: %f\n", prf);
  printf("prf_frac: %f\n", prf_frac);
  printf("nb: %f\n", nb);
  printf("nbg: %f\n", nbg);
  printf("nboff: %f\n", nboff);
  printf("bsl: %f\n", bsl);

  printf("kacoeff: %f, %f, %f\n", kacoeff[0], kacoeff[1], kacoeff[2]);
  printf("dopcoeff1: %f, %f, %f, %f\n", dopcoeff1[0], dopcoeff1[1], dopcoeff1[2], dopcoeff1[3]);
  printf("dopcoeff2: %f, %f, %f, %f\n", dopcoeff2[0], dopcoeff2[1], dopcoeff2[2], dopcoeff2[3]);


  if(nfilter % 2 != 1){
    fprintf(stderr, "filter length must be odd!\n");
    exit(1);
  }

  naz  = file_length(infp, nrg, sizeof(fcomplex));
  printf("file width: %d, file length: %d\n\n", nrg, naz);


  ka = array1d_float(nrg);
  dop1 = array1d_float(nrg);
  dop2 = array1d_float(nrg);

  nfa = array1d_float(nrg);
  freqs = array1d_float(nrg);
  freqe = array1d_float(nrg);
  bis = array1d_float(nrg);
  bie = array1d_float(nrg);
  bic = array1d_float(nrg);
  bica = array1d_float(nrg);

  in  = array2d_fcomplex(naz, nrg);
  out = array2d_fcomplex(naz, nrg);


  pri = 1.0/prf;
  nbc = nb + nbg;
  hnfilter = (nfilter - 1) / 2;

  //find burst starting line closest to first line and after first line
  for(i = -100000; i < 100000; i++){
    tmp1 = bsl + (nb + nbg) * i;
    if(tmp1 >= 0){
      bsl = tmp1;
      break;
    }
  }


  //calculate something
  for(i = 0; i < nrg; i++){
    
    //azimuth FM rate. we follow the convention ka > 0
    ka[i] = kacoeff[2] * i * i + kacoeff[1] * i + kacoeff[0];
    ka[i] = -ka[i];

    //doppler centroid frequency
    dop1[i] = dopcoeff1[0] + dopcoeff1[1] * i + dopcoeff1[2] * i * i + dopcoeff1[3] * i * i * i;
    dop1[i] *= prf;
    dop2[i] = dopcoeff2[0] + dopcoeff2[1] * i + dopcoeff2[2] * i * i + dopcoeff2[3] * i * i * i;
    dop2[i] *= prf;

    //full aperture length
    nfa[i] = prf * prf_frac / ka[i] / pri;

    //consider burst synchronization
    //these are the same for all columns
    if(fabs(nboff) >= 0.8 * nb){
      fprintf(stderr, "burst synchronization is too small!\n\n");
      exit(1);
    }
    if(nboff < 0){
      bsl_new = bsl - nboff;
    }
    else{
      bsl_new = bsl;
    }
    nb_new = nb - fabs(nboff);
    nbg_new = nbg + fabs(nboff);
    nbc_new = nbc;
    nbc_new_int = (int)(nbc_new + 0.5);

    //starting and ending doppler centroid frequency of the burst
    //if the overall doppler centroid frequency = 0
    freqs[i] = -(prf * prf_frac - nb_new * pri * ka[i]) / 2.0;
    freqe[i] =  (prf * prf_frac - nb_new * pri * ka[i]) / 2.0;

    //consider doppler centroid frequency
    freqs[i] += dop1[i];
    freqe[i] += dop1[i];

    //consider doppler centroid frequency of the other image
    tmp1 = dop2[i] - dop1[i];
    if(tmp1 > 0){
      freqs[i] += tmp1;
    }
    else{
      freqe[i] += tmp1;
    }

    //check if doppler centroid frequency difference too big
    if(freqe[i] - freqs[i] < nbc_new * pri * ka[i]){
      fprintf(stderr, "Doppler centroid frequency difference too large!\n\n");
      exit(1);
    }

    //starting and ending index of imaged area by the burst
    bic[i] = bsl_new + (nb_new - 1.0) / 2.0; //this should be the same for all columns
    bis[i] = freqs[i] / ka[i] / pri + bic[i];
    bie[i] = freqe[i] / ka[i] / pri + bic[i];
    bica[i] = (bis[i] + bie[i]) / 2.0;

  }


  //find the max and min of starting and ending index
  bis_min = bis[0];
  bis_max = bis[0];
  bie_min = bie[0];
  bie_max = bie[0];
  for(i = 0; i < nrg; i++){
    if(bis[i] < bis_min){
      bis_min = bis[i];
    }
    if(bis[i] > bis_max){
      bis_max = bis[i];
    }

    if(bie[i] < bie_min){
      bie_min = bie[i];
    }
    if(bie[i] > bie_max){
      bie_max = bie[i];
    }
  }

  
  //read in data
  readdata((fcomplex *)in[0], (size_t)naz * (size_t)nrg * sizeof(fcomplex), infp);

  //initialize output data
  for(j = 0; j < naz; j++){
    for(k = 0; k < nrg; k++){
      out[j][k].re = 0.0;
      out[j][k].im = 0.0;
    }
  }


  for(i = 0; i < nrg; i++){

    if((i + 1) % 100 == 0 || (i+1) == nrg)
      fprintf(stderr,"processing: %6d of %6d\r", i+1, nrg);
    if((i+1) == nrg)
      fprintf(stderr,"\n");

    
    //initialize start and ending line number
    if(dop1[i] > dop2[i]){
      bis_out = roundfi(bie[i]) + 1;
      //bie_out = roundfi(bie[i]) + 1 + (nbc_new - 1);

      //changed to use nbc_new_int. 27-JAN-2015
      bie_out = roundfi(bie[i]) + 1 + (nbc_new_int - 1);
    }
    else{
      bis_out = roundfi(bis[i]);
      //bie_out = roundfi(bis[i]) + (nbc_new - 1);

      //changed to use nbc_new_int. 27-JAN-2015
      bie_out = roundfi(bis[i]) + (nbc_new_int - 1);
    }

    //consider the filter length
    bis_in = bis_out - (nfilter - 1) / 2;
    bie_in = bie_out + (nfilter - 1) / 2;

    //to make circular convolution equivalent to linear convolution
    nfft = next_pow2(bie_in - bis_in + 1 + nfilter - 1);

    //initialize filter
    filter = array1d_fcomplex(nfft);
    filter_j = array1d_fcomplex(nfft);
    for(j = 0; j < nfft; j++){
      filter[j].re = 0.0;
      filter[j].im = 0.0;
    }


    nburst_new = (int)ceil(   fabs(freqe[i]-freqs[i]) / (nbc_new * pri * ka[i])   );

    //choose deramp center
    if(dop1[i] > dop2[i]){
      if(deramp_center_flag == 0){
        deramp_center = bic[i];
      }
      else if(deramp_center_flag == 1){
        deramp_center = (bica[i] + nbc_new);
      }
      else{
        deramp_center = bic[i] + (int)((nburst_new+1) / 2) * nbc_new;
      }
    }
    else{
      if(deramp_center_flag == 0){
        deramp_center = bic[i];
      }
      else if(deramp_center_flag == 1){
        deramp_center = bica[i];
      }
      else{
        deramp_center = bic[i] + (int)(nburst_new / 2) * nbc_new;
      }
    }

    //create filters
    for(j = 0; j <= nburst_new; j++){
      //center frequency of bandpass filter
      //determined by distance of raw burst center and deramp center
      if(dop1[i] > dop2[i]){
        bfc = (deramp_center - (bic[i] + j*nbc_new)) * pri * ka[i];
        
        //do not include first burst in this case
        if(j == 0){
          continue;
        }
      }
      else{
        bfc = (deramp_center - (bic[i] - j*nbc_new)) * pri * ka[i];
      
        //do not include last burst in this case
        if(j == nburst_new){
          break;
        }
      }

      //bandwidth of bandpass filter
      bfw = nb_new * pri * ka[i];

      //create filter: first sample corresponding to first fully convolution sample
      bandpass_filter(bfw/prf, bfc/prf, nfilter, nfft, nfilter-1, beta, filter_j);

      //add the filters to form the filter to be used
      for(k = 0; k < nfft; k++){
        filter[k].re += filter_j[k].re;
        filter[k].im += filter_j[k].im;
      }
    }

    //forward fft
    four1((float *)filter - 1, nfft, -1);

    //create deramp signal: this applies no matter whether dop1[i] is larger,
    //and no matter bic is on the left or right.
    deramp = array1d_fcomplex(nfft);
    for(j = 0; j < nfft; j++){
      //distance between fft center and deramp center
      //tmp1 = bis_in + (nfft - 1.0) / 2.0 - bic[i];
      tmp1 = bis_in + (nfft - 1.0) / 2.0 - deramp_center;

      //if(tmp1 <= 0){
      //  fprintf(stderr, "WARNING: very large doppler centroid frequnecy\n\n");
      //}
      //index used in deramp signal
      tmp2 = j - (nfft - 1.0) / 2.0 + tmp1;
      //deramp signal
      tmp3 = - PI * ka[i] * (tmp2 * pri) * (tmp2 * pri);
      deramp[j].re = cos(tmp3);
      deramp[j].im = sin(tmp3);
    }

    //rereamp signal
    reramp = array1d_fcomplex(nfft);
    for(j = 0; j < nfft; j++){
      reramp[j].re =  deramp[j].re;
      reramp[j].im = -deramp[j].im;
    }
    //circ_shift(reramp, nfft, -abs(nfilter-1));
    circ_shift(reramp, nfft, -abs(     (nfilter-1)/2     ));


/**********************************************/
/*             do the filtering               */ 
/**********************************************/


    //filter the data
    data = array1d_fcomplex(nfft);
    for(j = -10000; j < 10000; j++){
      //bis_out2 = bis_out + j * nbc_new;
      //bie_out2 = bie_out + j * nbc_new;
      //bis_in2  = bis_in  + j * nbc_new;
      //bie_in2  = bie_in  + j * nbc_new;
      
      //changed to use nbc_new_int. 27-JAN-2015
      bis_out2 = bis_out + j * nbc_new_int;
      bie_out2 = bie_out + j * nbc_new_int;
      bis_in2  = bis_in  + j * nbc_new_int;
      bie_in2  = bie_in  + j * nbc_new_int;

      //find data to be filtered
      if(bie_in2 <= -1){
        continue;
      }
      else if(bis_in2 >= naz){
        break;
      }
      else{
        //first zero the data
        for(k = 0; k < nfft; k++){
          data[k].re = 0.0;
          data[k].im = 0.0;
        }
        //get data
        for(k = bis_in2; k <= bie_in2; k++){
          if(k <= -1 || k >= naz){
            data[k-bis_in2].re = 0.0;
            data[k-bis_in2].im = 0.0;
          }
          else{
            data[k-bis_in2].re = in[k][i].re / sc;
            data[k-bis_in2].im = in[k][i].im / sc;
          }
        }
      }

      //deramp the data
      for(k = 0; k < nfft; k++){
        data[k] = cmul(data[k], deramp[k]);
      }

      //forward fft
      four1((float *)data - 1, nfft, -1);

      //multiplication in the frequency domain
      for(k = 0; k < nfft; k++)
        data[k] = cmul(data[k], filter[k]);

      //backward fft
      four1((float *)data - 1, nfft, 1);

      //reramp
      for(k = 0; k < nfft; k++){
        data[k] = cmul(data[k], reramp[k]);
      }

      //get the filtered data
      for(k = bis_out2; k <= bie_out2; k++){

        if(edgl_flag == 0){ //do not remove data on the edges
          edgl = 0;
        }
        else if(edgl_flag == 1){ //remove data less than half convolution
          edgl = (nfft - 1) / 2;
        }
        else{ //remove data of incomplete convolution
          edgl = nfft - 1;
        }

        if((k >= (0+edgl)) && (k <= naz-1-edgl)){
          out[k][i].re = data[k-bis_out2].re * sc / nfft;
          out[k][i].im = data[k-bis_out2].im * sc / nfft;
        }
      }

    }//j: block of data of each column
    free_array1d_fcomplex(filter);
    free_array1d_fcomplex(filter_j);
    free_array1d_fcomplex(deramp);
    free_array1d_fcomplex(reramp);
    free_array1d_fcomplex(data);
  }//i: each column

  writedata((fcomplex *)out[0], (size_t)naz * (size_t)nrg * sizeof(fcomplex), outfp);

  //free arrays
  free_array1d_float(ka);
  free_array1d_float(dop1);
  free_array1d_float(dop2);

  free_array1d_float(nfa);
  free_array1d_float(freqs);
  free_array1d_float(freqe);
  free_array1d_float(bis);
  free_array1d_float(bie);
  free_array1d_float(bic);
  free_array1d_float(bica);

  free_array2d_fcomplex(in);
  free_array2d_fcomplex(out);

  //close files
  fclose(infp);
  fclose(outfp);

  return 0;
}//end main()


//////////////////////////////////////////////////////////////////////////////////////////////////////



