#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dims1, *dims2, *cost_volume_dims;
    int dimx1, dimy1, dimx2, dimy2, numdims2, ch1, ch2;
    double * output, *min, *cost;
    int i,j,h,w,d;
  
    if (nrhs != 2)
    {
        mexErrMsgTxt("Two inputs required.");
    }
    
    //first is cost_volume
    cost_volume_dims = mxGetDimensions(prhs[0]);
    mwSize s1, s2, s3;
    s1 = cost_volume_dims[0];
    s2 = cost_volume_dims[1];
    s3 = cost_volume_dims[2];
    // second is minimum
    dims2 = mxGetDimensions(prhs[1]);
    dimy2 = (int)dims2[0]; dimx2 = (int)dims2[1];
    
    plhs[0] = mxCreateDoubleMatrix(s1, s2, mxREAL);
    cost = mxGetPr(prhs[0]);
    min = mxGetPr(prhs[1]);
    output = mxGetPr(plhs[0]);
    
    
    mexPrintf("s1=[%d], s2=[%d], s3=[%d], x_min=[%d], y_min=[%d]\n", s1, s2, s3, dimx2, dimy2);
    
    double sum = 0.0;
    double delta_d = 0.0;
    double delta_c = 0.0;
    double mean_c = 0;
    int ind_min = 0;
    double f1 = s3 / 3.0;
    double cost_min = 0;
    //do something
    double maxminnum =0.0;
    double maxdenom =0.0;
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    double tmp3 = 0.0; 
    double tmp4 = 0.0;
    
    for(i=0;i<s1;i++) //(i=hws+1;i<s1-hws-1;i++)
    {
        for(j=0;j<s2;j++) //(j=hws+1;j<s2-hws-1;j++)
        {
            //mexPrintf("pixel %d, %d\n", i,j);
            sum = 1.0;
            mean_c = 1.0;
            
            for (d = 0; d < s3; d++)
            {
                mean_c += cost[i + j * s1 + d * (s1 * s2)];
            }
            mean_c = mean_c / s3;
            //if (isnan(mean_c))
            //    mean_c = 0;
            //mexPrintf("mean_c=%.f, delta_c=%.3f\n",mean_c, delta_c);
            
            for (d = 0; d < s3; d++)
            {
                //mexPrintf("[%d,%d,%d]",i,j,d);
                ind_min = min[i + j * s1] - 1;
                cost_min = cost[i + j * s1 + ind_min * (s1 * s2)];
                delta_d = std::abs(d - ind_min);
                delta_c = cost[i + j * s1 + d * (s1 * s2)] - cost_min;
                tmp1 = cost[i + j * s1 + d * (s1 * s2)];
                tmp2 = ind_min;
                tmp3 = cost_min;
                tmp4 = mean_c/5.0;
                //mexPrintf("tmp1=%.3f, tmp2=%.3f, tmp3=%.3f, tmp4=%.3f\n",tmp1, tmp2, tmp3, tmp4);
                maxminnum = (pow(tmp3,2));
                maxdenom = (std::max(tmp4,1.0));
                
                //mexPrintf("maxminnum = %.3f, maxdenom = %.3f, div = %.3f, sum=%.3f\n",maxminnum, maxdenom, maxminnum/maxdenom, sum); //sum += delta_c;                
                sum += maxminnum / maxdenom;
            }
               
             
            //mexPrintf("sigma_cur[%f], sigma_aif[%f], sigma_cross[%f] \n",sigma_cur, sigma_aif, sigma_cross);
            output[i + j * s1] = i; // 1 / sum; 
            //mexPrintf("output %.3f", output[i + j * s1]);//-cur_fs[i*dimy1+j]; 
        }
    }
    

    return;
}