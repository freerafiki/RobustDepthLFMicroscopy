#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>
#include <vector>

int CalcMHWScore(std::vector<int> scores)
{
  size_t size = scores.size();

  if (size == 0)
  {
    return 0;  // Undefined, really.
  }
  else if (size == 1)
  {
    return scores[0];
  }
  else
  {
    sort(scores.begin(), scores.end());
    if (size % 2 == 0)
    {
      return (scores[size / 2 - 1] + scores[size / 2]) / 2;
    }
    else 
    {
      return scores[size / 2];
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dims1, *dims_sp_big, *dims_sp_small, *mask_dim;
    int d1, d2, d3, spb_y, spb_x, sps_y, sps_x, dimy4, dimx4;
    double * output_sp1, *output_sp2, *output_points, *cost_volume, *sp_big, *sp_small, *mask;
    int i,j,d,k;
  
    if (nrhs != 4)
    {
        mexErrMsgTxt("Four inputs required. \n1) Cost Volume \n2) Big Superpixels \n3) Smaller Superpixels \n4) Mask \n");
    }
    
    // first input is cost volume
    cost_volume = mxGetPr(prhs[0]);
    dims1 = mxGetDimensions(prhs[0]);
    d1 = (int)dims1[0]; d2 = (int)dims1[1]; d3 = (int)dims1[2];
    // second and third are superpixels 
    sp_big = mxGetPr(prhs[1]);
    dims_sp_big = mxGetDimensions(prhs[1]);
    spb_y = (int)dims_sp_big[0]; spb_x = (int)dims_sp_big[1];
    sp_small = mxGetPr(prhs[2]);
    dims_sp_small = mxGetDimensions(prhs[2]);
    sps_y = (int)dims_sp_small[0]; sps_x = (int)dims_sp_small[1];
    // fourth is the mask
    mask = mxGetPr(prhs[3]);
    mask_dim = mxGetDimensions(prhs[3]);
    dimy4 = (int)mask_dim[0]; dimx4 = (int)mask_dim[1];

    // maximum superpixel 1
    const int max_sp_big = sp_big[d2*(d2-1)+d1-1];
    // maximum superpixel 2
    const int max_sp_small = sp_small[d2*(d2-1)+d1-1];
    
    // output will be three matrices with the argmin (depth)
    plhs[0] = mxCreateDoubleMatrix(d1, d2, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(d1, d2, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(d1, d2, mxREAL);
    output_sp1 = mxGetPr(plhs[0]);
    output_sp2 = mxGetPr(plhs[1]);
    output_points = mxGetPr(plhs[2]);
    
    // to fill the superpixels parts
    std::vector< std::vector<int> > sp_big_values;
    std::vector< std::vector<int> > sp_small_values;
    for (k=0;k<max_sp_big;k++)
        sp_big_values.push_back(std::vector<int>());
    for (k=0;k<max_sp_small;k++)
        sp_small_values.push_back(std::vector<int>());
    int sp_big_median [max_sp_big];
    int sp_small_median [max_sp_small];
    
    double min_value = 10000.0;
    int min_index = 0;
    
    //loop on the image
    for(i=0;i<d1;i++) // 300; i < 350; i++) //
    {
        for(j=0;j<d2;j++) //300; j < 350; j++) // 
        {
            if (mask[i*d1+j] > 0)
            {
                //mexPrintf("i=%d, j=%d\n", i,j);
                min_value = 100000.0;
                min_index = 0;
                for (d=0;d<d3;d++)
                {
                    if (cost_volume[i + j * d1 + d * (d1 * d2)] < min_value)
                    {
                        min_value = cost_volume[i + j * d1 + d * (d1 * d2)];
                        min_index = d;
                    }
                }
                sp_big_values[(int)(sp_big[i*d2+j])].push_back(min_index);
                sp_small_values[(int)(sp_small[i*d2+j])].push_back(min_index);
                output_points[i*d2+j] = min_index; //-cur_fs[i*dimy1+j];
            }
            else
            {
                output_points[i*d2+j] = -1;
            }
        }
    }
    //mexPrintf("second loop");
    // loop over arrays to create median 
    for (k=0;k<max_sp_big;k++)
    {
        sp_big_median[k] = CalcMHWScore((sp_big_values[k]));
    }
    for (k=0;k<max_sp_small;k++)
    {
        sp_small_median[k] = CalcMHWScore((sp_small_values[k]));
    }
    //fill the sp map
    for(i=0;i<d1;i++) // 300; i < 350; i++) //
    {
        for(j=0;j<d2;j++) //300; j < 350; j++) // 
        {
            if (mask[i*d2+j] > 0)
            {
                output_sp1[i*d2+j] = sp_big_median[(int)(sp_big[i*d2+j])];
                output_sp2[i*d2+j] = sp_small_median[(int)(sp_small[i*d2+j])];
            }
            else
            {
                output_sp1[i*d2+j] = -1;
                output_sp2[i*d2+j] = -1;
            }
        }
    }
    return;
}
