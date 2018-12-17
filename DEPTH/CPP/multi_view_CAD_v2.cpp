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
    const mwSize *stack_dims;
    int d1, d2, d3;
    double *sad, *census, *big_mat, *cost_volume;
    int i,j,h,w,d,dd,ci, dmin, dmax, disp, ws, hws;
    double alpha, census_hd, sum;
  
    if (nrhs != 5)
    {
        mexErrMsgTxt("Five inputs required. Stack of images, minimum disparity, maximum disparity, window size and alpha.");
    }
    
    //figure out dimensions
    stack_dims = mxGetDimensions(prhs[0]);
    d1 = (int)stack_dims[0]; // y dimension
    d2 = (int)stack_dims[1]; // x dimension
    d3 = (int)stack_dims[2]; // number of images
    dmin = *mxGetPr(prhs[1]);
    dmax = *mxGetPr(prhs[2]);
    disp = (dmax - dmin)+1;
    ws = *mxGetPr(prhs[3]); // window size
    hws = floor(ws/2.0);
    alpha = *mxGetPr(prhs[4]); // alpha to combine census and sad
    // build a cost volume that has 4 dimensions: x,y, num of disparities and num of images
    // it contains at every image the cost of the matching central image against i-th iamge
    // in the central image the cost of all images summed (multi-view cost)
    const mwSize cost_volume_dims[]={300,300,disp,d3}; //d1,d2,disp,d3};
    plhs[0] = mxCreateNumericArray(4, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    big_mat = mxGetPr(prhs[0]);
    cost_volume = mxGetPr(plhs[0]);
    mwSize s1, s2, s3, s4;
    s1 = 300; //cost_volume_dims[0]; //x 
    s2 = 300; //cost_volume_dims[1]; //y
    s3 = disp; //cost_volume_dims[2]; //disparities
    s4 = d3; //cost_volume_dims[3]; //num of images
    // variables for computing costs
    int c_im = floor(d3/2.0); // central image index
    sum = 0.0;
    census_hd = 0.0;
    double current_img_sum = 0.0;
    double current_img_census = 0.0;
    double TRUNC_AD = d3*50;
    int left_census[ws][ws];
    int right_census[ws][ws];
    int c1=0; 
    int c2=0; 
    int i_=0, j_=0;
    int pad = std::max(std::abs(dmax), std::abs(dmin));
    mexPrintf("dmin=%d, dmax=%d, ws=%d, pad=%d, alpha=%1.1f \n", dmin, dmax, ws, pad, alpha);
    mexPrintf("d1=%d, d2=%d, d3=%d, d4=%d\n", d1, d2, disp, d3);
    mexPrintf("d1=%d, d2=%d, d3=%d, d4=%d\n", s1, s2, s3, s4);
    //do something
    for(i=pad+hws+1;i<s1-pad-hws-1;i++) //320; i<370; i++) //(i=hws+1;i<s1-hws-1;i++)
    {
        //if ((i%5) == 0)
        //    mexPrintf("ind=%d\n",i+(j)*s1+ci*(s1*s2));
        for(j=pad+hws+1;j<s2-pad-hws-1;j++) //320; j<370; j++) //(j=hws+1;j<s2-hws-1;j++)
        {
            sum = 0.0;
            census_hd = 0.0;
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    if (big_mat[i+h+(j+w)*s1+ci*(s1*s2)] < big_mat[i+(j)*s1+ci*(s1*s2)])
                        left_census[h+hws][w+hws] = 0;
                    else
                        left_census[h+hws][w+hws] = 1;
                }
            }
            for (d = dmin; d <= dmax; d++)
            {
                dd = d - dmin;
                sum = 0.0;
                census_hd = 1;
                if (d > -0.1 && d < 0.1)
                {
                    sum = 1;
                    census_hd = 1;
                }
                else
                {
                    sum = 0.0;
                    census_hd = 0.0;
                    for (ci = 0; ci < d3; ci++)
                    {
                        map_index_to_coords(ci, d3, c1, c2);
                        i_ = round(i+c1*0.866*(d));
                        j_ = round(j+c2*0.5*(d));

                        current_img_census = 0.0;
                        current_img_sum = 0.0;
                        for(h=-hws;h<hws;h++)
                        {
                            for(w=-hws;w<hws;w++)
                            {
                                if (big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] < big_mat[i_+(j_)*s1+ci*(s1*s2)])
                                    right_census[h+hws][w+hws] = 0;
                                else
                                    right_census[h+hws][w+hws] = 1;
                                current_img_census += std::abs(right_census[h+hws][w+hws] - left_census[h+hws][w+hws]);
                                current_img_sum += std::min(std::abs(big_mat[i_+h+(j_+w)*s1+ci*(s1*s2)] - big_mat[i+h+(j+w)*s1+c_im*(s1*s2)]), TRUNC_AD); 
                            }
                        }

                        //sum /= 1;
                        current_img_sum /= (ws*ws*TRUNC_AD);
                        current_img_census /= (ws*ws);
                        sum += current_img_sum;
                        census_hd += current_img_census;
                        // cost_volume[:,:,dd,i] = is the cost volume for disparity d = (dd - dmin) for 
                        // central image against i-th image
                        //mexPrintf("i = %d, j=%d, dd=%d, ci=%d, ind=%d\n",i,j,dd,ci,i + j * s1 + dd * (s1 * s2) + ci * (s1 * s2 * s3));
                        cost_volume[i + j * s1 + dd * (s1 * s2) + ci * (s1 * s2 * s3)] = alpha*current_img_sum + (1-alpha)*current_img_census;
                        //mexPrintf("ci=%d, c1=%d, c2=%d, i=%d, j=%d, i_=%d, j_=%d, d=%d\n",ci, c1, c2, i, j, i_, j_, d);
                    }
                }
                
                //mexPrintf("i = %d, j=%d, dd=%d, c_im=%d, ind=%d\n",i,j,dd,c_im,i + j * s1 + dd * (s1 * s2) + c_im * (s1 * s2 * s3));
                cost_volume[i + j * s1 + dd * (s1 * s2) + c_im * (s1 * s2 * s3)] = alpha*sum + (1-alpha)*census_hd;
            }
        }
    }
    return;
}


