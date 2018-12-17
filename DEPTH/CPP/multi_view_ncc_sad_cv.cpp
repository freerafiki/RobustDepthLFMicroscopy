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
        mexErrMsgTxt("Five inputs required. Central image, matrices with other images, disparities vector, window size and alpha.");
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
    double alpha = *mxGetPr(prhs[4]);
    int ws = *mxGetPr(prhs[3]);
    int hws = floor(ws/2.0);
    int c_im = floor(s3/2.0);
    double ncc_sum = 0.0;
    double sad_sum = 0.0;
    int c1=0; 
    int c2=0; 
    int i_=0, j_=0;
    int half_disp = floor(dimx3/2);
    int dmin = 25;
    int dmax = 35;
    double mean_cen_img = 0.0;
    double mean_cur_img = 0.0;
    double tmp_cen_img = 0.0;
    double tmp_cur_img = 0.0;
    double sigma_cur_img = 0.0;
    double sigma_cen_img = 0.0;
    double sigma_cross = 0.0;
    double disp_sad = 0.0;
    double disp_ncc = 0.0;
    //mexPrintf("dmin=%d, dmax=%d \n", dmin, dmax);
    //do something
    for(i=dmax+hws+1;i<s1-dmax-hws-1;i++) // 100;i<200;i++) // //(i=hws+1;i<s1-hws-1;i++)
    {
        for(j=dmax+hws+1;j<s2-dmax-hws-1;j++) //  300;j<400;j++) //  (j=hws+1;j<s2-hws-1;j++)
        {
            // calculate mean for the reference (central) image
            mean_cen_img = 0.0;
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    mean_cen_img += c_img[i+h+(j+w)*s1]; 
                }
            } 
            mean_cen_img /= (ws*ws);
            //if (mean_cen_img > 0)
            //    mexPrintf("mci=%.3f\n",mean_cen_img);
            for (d = 0; d < s3; d++)
            {
                // these contains diff values for disparity d
                disp_sad = 0.0;
                disp_ncc = 0.0;
                disp = d - half_disp;
                for (ci = 0; ci < ch2; ci++)
                {
                    // these values are for one image
                    ncc_sum = 0.0;
                    sad_sum = 0.0;
                    // get coords for the each image
                    // for each ci we are using one of the EIs
                    map_index_to_coords(ci, ch2, c1, c2);
                    i_ = round(i+c1*0.866*(d));
                    j_ = round(j+c2*0.5*(d));
                    
                    // calculate mean 
                    mean_cur_img = 0;
                    for(h=-hws;h<hws;h++)
                    {
                        for(w=-hws;w<hws;w++)
                        {
                            mean_cur_img += big_mat[i+h+(j+w)*s1+ci*(s1*s2)]; 
                        }
                    } 
                    mean_cur_img /= (ws*ws);
                    
                    // calculate NCC and SAD for the window
                    sigma_cur_img = 0.0;
                    sigma_cen_img = 0.0;
                    sigma_cross = 0.0;
                    for(h=-hws;h<hws;h++)
                    {
                        for(w=-hws;w<hws;w++)
                        {
                            sad_sum += pow(big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - big_mat[i+h+(j+w)*s1+4*(s1*s2)],2);
                            tmp_cur_img = (big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - mean_cur_img);
                            sigma_cur_img += pow(tmp_cur_img,2);
                            tmp_cen_img = (c_img[i+h+(j+w)*s1] - mean_cen_img);
                            sigma_cen_img += pow(tmp_cen_img,2);
                            sigma_cross += (tmp_cur_img)*(tmp_cen_img);
                            /*
                            sad_sum += (big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - big_mat[i+h+(j+w)*s1+c_im*(s1*s2)])*(big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - big_mat[i+h+(j+w)*s1+c_im*(s1*s2)]); 
                            tmp_cur_img = (pow(big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - mean_cur_img,2));
                            sigma_cur_img += tmp_cur_img;
                            tmp_cen_img = (pow(c_img[i+h+(j+w)*s1] - mean_cen_img,2));
                            sigma_cen_img += tmp_cen_img;
                            sigma_cross += (tmp_cur_img)*(tmp_cen_img);
                             */
                        }
                    }
                    sigma_cen_img /= (ws*ws);
                    sigma_cur_img /= (ws*ws);
                    sigma_cross /= (ws*ws);
                    //normalize
                    //sad_sum /= (ws*ws);
                    sad_sum /= (ws*ws); //*(big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - big_mat[i+h+(j+w)*s1+c_im*(s1*s2)]); 
                    ncc_sum = (sigma_cur_img*sigma_cen_img) / (sigma_cross);
                    //mexPrintf("ci=%d, c1=%d, c2=%d, i=%d, j=%d, i_=%d, j_=%d, d=%d\n",ci, c1, c2, i, j, i_, j_, d);
                    disp_sad += sad_sum;
                    disp_ncc += ncc_sum;
                }
                output[i + j * s1 + d * (s1 * s2)] = alpha * disp_ncc + (1.0-alpha) * disp_sad; 
            }
        }
    }
    return;
}


