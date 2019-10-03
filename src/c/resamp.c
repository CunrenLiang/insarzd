//////////////////////////////////////
// Cunren Liang, NASA JPL/Caltech
// Copyright 2017
//////////////////////////////////////


#include "resamp.h"

int main(int argc, char *argv[]){
  FILE *slc2fp;
  FILE *rslc2fp;
  FILE *rgoffp;
  FILE *azoffp;
  int rgflag;
  int azflag;
  int nrg1;
  int naz1;
  int nrg2;
  int naz2;
  float prf;
  float dopcoeff[4];
  float rgcoef[10];
  float azcoef[10];
  float azpos_off;
  float beta;
  int n;
  int m;
  int interp_method;
  int edge_method;
  float rgpos;
  float azpos;
  float rgoff;
  float azoff;
  float rgoff1;
  float azoff1;
  float *rgoff2;
  float *azoff2;
  float rg2;
  float az2;
  int rgi2;
  int azi2;
  float rgf;
  float azf;
  int rgfn;
  int azfn;
  int hnm;
  int hn;
  float *sincc;
  float *kaiserc;
  float *kernel;
  float *rgkernel;
  float *azkernel;
  fcomplex *azkernel_fc;
  fcomplex *rgrs;
  fcomplex *azca;
  fcomplex *rgrsb;
  fcomplex *azrs;
  float *dop;
  float dopx;
  fcomplex **inb;
  int i, j, k, k1, k2;
  int tmp1, tmp2;
  int zero_flag;
  float ftmp1, ftmp2;
  fcomplex fctmp1, fctmp2;
  beta = 2.5;
  n = 9;
  m = 10000;
  interp_method = 0;
  edge_method = 0;
  if(argc < 14){
    fprintf(stderr, "\nusage: %s, slc2 rslc2 rgoff azoff nrg1 naz1 nrg2 naz2 prf dopcoeff[0]-[3] <rgcoef[0]-[9]> <azcoef[0]-[9]> <azpos_off>", argv[0]);
    fprintf(stderr, "\nmandatory:\n");
    fprintf(stderr, "slc2:  slave image\n");
    fprintf(stderr, "rslc2: resampled slave image\n");
    fprintf(stderr, "rgoff: range offset file. if no range offset file, specify fake\n");
    fprintf(stderr, "azoff: azimuth offset file. if no azimuth offset file, specify fake\n");
    fprintf(stderr, "nrg1:  number of columns in master image\n");
    fprintf(stderr, "naz1:  number of lines in master image\n");
    fprintf(stderr, "nrg2:  number of columns in slave image\n");
    fprintf(stderr, "naz2:  number of lines in slave image\n");
    fprintf(stderr, "prf:   PRF of slave image\n");
    fprintf(stderr, "dopcoeff[0]-[3]: Doppler centroid frequency coefficents\n");
    fprintf(stderr, "\noptional:\n");
    fprintf(stderr, "rgcoef[0]-[9]:   range offset polynomial coefficents. First of two fit results of resamp_roi\n");
    fprintf(stderr, "azcoef[0]-[9]:   azimuth offset polynomial coefficents. First of two fit results of resamp_roi\n");
    fprintf(stderr, "azpos_off:       azimuth position offset. Azimuth line number (column 3) of first offset in culled offset file\n\n");
    exit(1);
  }
  slc2fp = openfile(argv[1], "rb");
  rslc2fp = openfile(argv[2], "wb");
  rgflag = 0;
  azflag = 0;
  if (strcmp(argv[3], "fake") == 0){
    rgflag = 0;
    printf("range offset file not provided\n");
  }
  else{
    rgflag = 1;
    rgoffp = openfile(argv[3], "rb");
  }
  if (strcmp(argv[4], "fake") == 0){
    azflag = 0;
    printf("azimuth offset file not provided\n");
  }
  else{
    azflag = 1;
    azoffp = openfile(argv[4], "rb");
  }
  nrg1 = atoi(argv[5]);
  naz1 = atoi(argv[6]);
  nrg2 = atoi(argv[7]);
  naz2 = atoi(argv[8]);
  prf = atof(argv[9]);
  for(i = 0; i < 4; i++){
    dopcoeff[i] = atof(argv[10+i]);
  }
  for(i = 0; i < 10; i++){
    if(argc > 14 + i)
      rgcoef[i] = atof(argv[14+i]);
    else
      rgcoef[i] = 0.0;
  }
  for(i = 0; i < 10; i++){
    if(argc > 24 + i)
      azcoef[i] = atof(argv[24+i]);
    else
      azcoef[i] = 0.0;
  }
  if(argc > 34)
    azpos_off = atof(argv[34]);
  printf("\n\ninput parameters:\n");
  printf("slc2: %s\n", argv[1]);
  printf("rslc2: %s\n", argv[2]);
  printf("rgoff: %s\n", argv[3]);
  printf("azoff: %s\n\n", argv[4]);
  printf("nrg1: %d\n", nrg1);
  printf("naz1: %d\n", naz1);
  printf("nrg2: %d\n", nrg2);
  printf("naz2: %d\n\n", naz2);
  printf("prf: %f\n\n", prf);
  for(i = 0; i < 4; i++){
    printf("dopcoeff[%d]: %e\n", i, dopcoeff[i]);
  }
  printf("\n");
  for(i = 0; i < 10; i++){
    printf("rgcoef[%d]: %e\n", i, rgcoef[i]);
  }
  printf("\n");
  for(i = 0; i < 10; i++){
    printf("azcoef[%d]: %e\n", i, azcoef[i]);
  }
  printf("\n");
  printf("azpos_off: %f\n\n", azpos_off);
hn = n / 2;
hnm = n * m / 2;
rgoff2 = array1d_float(nrg1);
azoff2 = array1d_float(nrg1);
sincc = vector_float(-hnm, hnm);
kaiserc = vector_float(-hnm, hnm);
kernel = vector_float(-hnm, hnm);
rgkernel = vector_float(-hn, hn);
azkernel = vector_float(-hn, hn);
azkernel_fc = vector_fcomplex(-hn, hn);
rgrs = vector_fcomplex(-hn, hn);
azca = vector_fcomplex(-hn, hn);
rgrsb = vector_fcomplex(-hn, hn);
azrs = array1d_fcomplex(nrg1);
dop = array1d_float(nrg2);
inb = array2d_fcomplex(naz2, nrg2);
sinc(n, m, sincc);
kaiser(n, m, kaiserc, beta);
for(i = -hnm; i <= hnm; i++)
  kernel[i] = kaiserc[i] * sincc[i];
printf("\n");
for(i = 0; i < nrg2; i++){
  dop[i] = dopcoeff[0] + dopcoeff[1] * i + dopcoeff[2] * i * i + dopcoeff[3] * i * i * i;
  dop[i] *= prf;
  if(i % 500 == 0)
    printf("range sample: %5d, doppler centroid frequency: %8.2f Hz\n", i, dop[i]);
}
printf("\n");
readdata((fcomplex *)inb[0], (size_t)naz2 * (size_t)nrg2 * sizeof(fcomplex), slc2fp);
for(i = 0; i < naz1; i++){
  if((i + 1) % 100 == 0)
    fprintf(stderr,"processing line: %6d of %6d\r", i+1, naz1);
    if (rgflag == 1){
      readdata((float *)rgoff2, nrg1 * sizeof(float), rgoffp);
    }
    if (azflag == 1){
      readdata((float *)azoff2, nrg1 * sizeof(float), azoffp);
    }
  for(j = 0; j < nrg1; j++){
    azrs[j].re = 0.0;
    azrs[j].im = 0.0;
  }
  for(j = 0; j < nrg1; j++){
    rgpos = j;
    azpos = i - azpos_off;
    rgoff1 = rgcoef[0] + azpos*(rgcoef[2] + \
             azpos*(rgcoef[5] + azpos*rgcoef[9])) + \
             rgpos*(rgcoef[1] + rgpos*(rgcoef[4] + \
             rgpos*rgcoef[8])) + \
             rgpos*azpos*(rgcoef[3] + azpos*rgcoef[6] + \
             rgpos*rgcoef[7]);
    azoff1 = azcoef[0] + azpos*(azcoef[2] + \
             azpos*(azcoef[5] + azpos*azcoef[9])) + \
             rgpos*(azcoef[1] + rgpos*(azcoef[4] + \
             rgpos*azcoef[8])) + \
             rgpos*azpos*(azcoef[3] + azpos*azcoef[6] + \
             rgpos*azcoef[7]);
    if (rgflag == 1){
      rgoff = rgoff1 + rgoff2[j];
    }
    else{
      rgoff = rgoff1;
    }
    if (azflag == 1){
      azoff = azoff1 + azoff2[j];
    }
    else{
      azoff = azoff1;
    }
    rg2 = j + rgoff;
    az2 = i + azoff;
    rgi2 = roundfi(rg2);
    azi2 = roundfi(az2);
    rgf = rg2 - rgi2;
    azf = az2 - azi2;
    rgfn = roundfi(rgf * m);
    azfn = roundfi(azf * m);
    for(k = -hn; k <= hn; k++){
      tmp1 = k * m - rgfn;
      tmp2 = k * m - azfn;
      if(tmp1 > hnm) tmp1 = hnm;
      if(tmp2 > hnm) tmp2 = hnm;
      if(tmp1 < -hnm) tmp1 = -hnm;
      if(tmp2 < -hnm) tmp2 = -hnm;
      rgkernel[k] = kernel[tmp1];
      azkernel[k] = kernel[tmp2];
    }
    for(k1 = -hn; k1 <= hn; k1++){
      rgrs[k1].re = 0.0;
      rgrs[k1].im = 0.0;
      if(edge_method == 0){
        if(azi2 < hn || azi2 > naz2 - 1 - hn || rgi2 < hn || rgi2 > nrg2 - 1 - hn){
          continue;
        }
      }
      else if(edge_method == 1){
        if(azi2 < 0 || azi2 > naz2 - 1 || rgi2 < 0 || rgi2 > nrg2 - 1){
          continue;
        }
      }
      else{
        if(azi2 < -hn || azi2 > naz2 - 1 + hn || rgi2 < -hn || rgi2 > nrg2 - 1 + hn){
          continue;
        }
      }
      for(k2 = -hn; k2 <= hn; k2++){
        if(azi2 + k1 < 0 || azi2 + k1 > naz2 - 1 || rgi2 + k2 < 0 || rgi2 + k2 > nrg2 - 1)
          continue;
        rgrs[k1].re += inb[azi2 + k1][rgi2 + k2].re * rgkernel[k2];
        rgrs[k1].im += inb[azi2 + k1][rgi2 + k2].im * rgkernel[k2];
      }
    }
    for(k = -hn; k <= hn; k++){
      if(rgrs[k].re == 0.0 && rgrs[k].im == 0.0)
        continue;
      dopx = dopcoeff[0] + dopcoeff[1] * rg2 + dopcoeff[2] * rg2 * rg2 + dopcoeff[3] * rg2 * rg2 * rg2;
      dopx *= prf;
      ftmp1 = 2.0 * PI * dopx * k / prf;
      azca[k].re = cos(ftmp1);
      azca[k].im = sin(ftmp1);
      if(interp_method == 0){
       rgrsb[k] = cmul(rgrs[k], cconj(azca[k]));
        azrs[j].re += rgrsb[k].re * azkernel[k];
        azrs[j].im += rgrsb[k].im * azkernel[k];
      }
      else{
        azkernel_fc[k].re = azca[k].re * azkernel[k];
        azkernel_fc[k].im = azca[k].im * azkernel[k];
        azrs[j] = cadd(azrs[j], cmul(rgrs[k], azkernel_fc[k]));
      }
    }
    if(interp_method == 0){
      ftmp1 = 2.0 * PI * dopx * azf / prf;
      fctmp1.re = cos(ftmp1);
      fctmp1.im = sin(ftmp1);
      azrs[j] = cmul(azrs[j], fctmp1);
    }
  }
  writedata((fcomplex *)azrs, nrg1 * sizeof(fcomplex), rslc2fp);
}
fprintf(stderr,"processing line: %6d of %6d\n", naz1, naz1);
free_array1d_float(rgoff2);
free_array1d_float(azoff2);
free_vector_float(sincc, -hnm, hnm);
free_vector_float(kaiserc, -hnm, hnm);
free_vector_float(kernel, -hnm, hnm);
free_vector_float(rgkernel, -hn, hn);
free_vector_float(azkernel, -hn, hn);
free_vector_fcomplex(azkernel_fc, -hn, hn);
free_vector_fcomplex(rgrs, -hn, hn);
free_vector_fcomplex(azca, -hn, hn);
free_vector_fcomplex(rgrsb, -hn, hn);
free_array1d_fcomplex(azrs);
free_array1d_float(dop);
free_array2d_fcomplex(inb);
fclose(slc2fp);
fclose(rslc2fp);
if (rgflag == 1){
  fclose(rgoffp);
}
if (azflag == 1){
  fclose(azoffp);
}
return 0;
}
