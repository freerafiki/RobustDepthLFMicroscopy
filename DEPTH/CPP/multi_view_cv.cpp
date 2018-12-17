#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void map_index_to_coords(int ci, int ch, int &c1, int &c2) {
  
  if (ci == 0)
  {
      c1 = -1;
      c2 = -1;
  }
  else if (ci == 1)
  {
      c1 = -1;
      c2 = +1;
  }
  else if (ci == 2)
  {
      c1 = 0;
      c2 = -2;
  }
  else if (ci == 3)
  {
      c1 = 0;
      c2 = 0;
  }
  else if (ci == 4)
  {
      c1 = 0;
      c2 = +2;
  }
  else if (ci == 5)
  {
      c1 = +1;
      c2 = -1;
  }
  else if (ci == 6)
  {
      c1 = +1;
      c2 = +1;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dims1, *dims2, *dims3;
    int dimx1, dimy1, numdims1, dimx2, dimy2, numdims2, ch1, ch2, dimx3, dimy3;
    double * output, *c_img, *big_mat, *disparities;
    int i,j,h,w,d, disp, ci;
    int num_disp;
  
    if (nrhs != 5)
    {
        mexErrMsgTxt("Five inputs required. Central image, matrices with other images, disparities vector, dmin and dmax.");
    }
    
    //figure out dimensions
    dims1 = mxGetDimensions(prhs[0]);
    numdims1 = mxGetNumberOfDimensions(prhs[0]);
    dimy1 = (int)dims1[0]; dimx1 = (int)dims1[1];
    if (numdims1 == 3)
        ch1 = (int)dims1[2];
    else
        ch1 = 1;
    dims2 = mxGetDimensions(prhs[1]);
    numdims2 = mxGetNumberOfDimensions(prhs[1]);
    dimy2 = (int)dims2[0]; dimx2 = (int)dims2[1];
    if (numdims2 == 3)
        ch2 = (int)dims2[2];
    else
        ch2 = 1;
    dims3 = mxGetDimensions(prhs[2]);
    dimy3 = (int)dims3[0]; dimx3 = (int)dims3[1];
    
    //mexPrintf("h1=[%d], w1=[%d], ch1=[%d], h2=[%d], w2=[%d], ch2=[%d], h3=[%d], w3=[%d]\n", dimx1, dimy1, ch1, dimx2, dimy2, ch2, dimx3, dimy3);
    
    const mwSize cost_volume_dims[]={dimx1,dimy1,dimx3};
    plhs[0] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    c_img = mxGetPr(prhs[0]);
    big_mat = mxGetPr(prhs[1]);
    disparities = mxGetPr(prhs[2]);
    output = mxGetPr(plhs[0]);
    mwSize s1, s2, s3;
    s1 = cost_volume_dims[0];
    s2 = cost_volume_dims[1];
    s3 = cost_volume_dims[2];
    int diameter = s1;
    int ws = 41;
    int hws = floor(ws/2.0);
    int mid = floor(s1/2.0);
    int c_im = floor(ch2/2.0);
    double sum = 0.0;
    int c1=0; 
    int c2=0; 
    int i_, j_;
    int half_disp = floor(dimx3/2);
    int dmin = *mxGetPr(prhs[3]);
    int dmax = *mxGetPr(prhs[4]);
    //mexPrintf("dmin=%d, dmax=%d \n", dmin, dmax);
    //do something
    for(i=dmax+hws+1;i<s1-dmax-hws-1;i++) //(i=hws+1;i<s1-hws-1;i++)
    {
        for(j=dmax+hws+1;j<s2-dmax-hws-1;j++) //(j=hws+1;j<s2-hws-1;j++)
        {
            sum = 0.0;
            for (d = 0; d < s3; d++)
            {
                sum = 0.0;
                disp = d - half_disp;
                for (ci = 0; ci < ch2; ci++)
                {
                    map_index_to_coords(ci, ch2, c1, c2);
                    i_ = round(i+c1*0.866*(d));
                    j_ = round(j+c2*0.5*(d));
                    
                    for(h=-hws;h<hws;h++)
                    {
                        for(w=-hws;w<hws;w++)
                        {
                            sum += big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - big_mat[i+h+(j+w)*s1+c_im*(s1*s2)]; 
                        }
                    }
                    
                    sum += big_mat[i_+(j_)*s1+ci*(s1*s2)] - big_mat[i+(j)*s1+c_im*(s1*s2)];
                    //mexPrintf("ci=%d, c1=%d, c2=%d, i=%d, j=%d, i_=%d, j_=%d, d=%d\n",ci, c1, c2, i, j, i_, j_, d);
                    
                // - big_mat[i + j * s1 + ci * (s1 * s2)])*(c_img[i + j * s1] - big_mat[i + j * s1 + ci * (s1 * s2)]);
                }
                output[i + j * s1 + d * (s1 * s2)] = sum; 
            }
        }
    }
    return;
}


