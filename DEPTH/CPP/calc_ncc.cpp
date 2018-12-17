#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dims1, *dims2;
    int dimx1, dimy1, numdims1, dimx2, dimy2, numdims2, ch1, ch2;
    double * output, *aif, *cur_fs;
    int i,j,h,w;
  
    if (nrhs != 3)
    {
        mexErrMsgTxt("Three inputs required.");
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
        ch2 = (int)dims1[2];
    else
        ch2 = 1;
    
    //mexPrintf("h1=[%d], w1=[%d], h2=[%d], w2=[%d]\n", dimx1, numdims1, dimx2, numdims2);
    
    plhs[0] = mxCreateDoubleMatrix(dimx1, dimy1, mxREAL);
    cur_fs = mxGetPr(prhs[0]);
    aif = mxGetPr(prhs[1]);
    output = mxGetPr(plhs[0]);
    int ws = *mxGetPr(prhs[2]);
    int hws = floor(ws/2.0);
    
    double mean_p1 = 0.0;
    double mean_p2 = 0.0;
    double mean_p12 = 0.0;
    double sigma_p1 = 0.0;
    double sigma_p2 = 0.0;
    double ncc_val = 0.0;
    double tmp_p1 = 0.0;
    double tmp_p2 = 0.0;
    //do something
    for(i=hws+1;i<dimx1-hws-1;i++) // 250; i<300; i++) //
    {
        for(j=hws+1;j<dimy1-hws-1;j++) // 250; j<300; j++) //
        {
            // mean of the window
            mean_p1 = 0.0;
            mean_p2 = 0.0;
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    mean_p1 += aif[(i+h)*dimy1+(j+w)];
                    mean_p2 += cur_fs[(i+h)*dimy1+(j+w)];
                }
            }
            mean_p1 /= (ws*ws);
            mean_p2 /= (ws*ws);
            
            // calcualte NCC
            sigma_p1 = 0.0;
            sigma_p2 = 0.0;
            mean_p12 = 0.0;
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    tmp_p1 = (aif[(i+h)*dimy1+(j+w)] - mean_p1);
                    tmp_p2 = (cur_fs[(i+h)*dimy1+(j+w)] - mean_p2);
                    mean_p12 += tmp_p1*tmp_p2;
                    sigma_p1 += pow(tmp_p1,2);
                    sigma_p2 += pow(tmp_p2,2);
                }
            }   
            sigma_p1 /= (ws*ws);
            sigma_p2 /= (ws*ws);
            mean_p12 /= (ws*ws);
            //mexPrintf("sigma_cur[%f], sigma_aif[%f], sigma_cross[%f] \n",sigma_p1, sigma_p1, (mean_p12) / (sigma_p1 * sigma_p2));
            output[i*dimy1+j] = (mean_p12) / (sigma_p1 * sigma_p2); //-cur_fs[i*dimy1+j]; 
        }
    }
    /* allocate the matrix to be returned */
    //plhs[0] = mxCreateDoubleMatrix(h2, w2, mxREAL);
    //diff = mxGetPr(plhs[0]);
    //for(int i=0;i<h2*w2;i++) 
    //    *(diff+i)=*(img1_ptr+i)*2.0;
    /*
     * int patch_size = 3;
    for(int k=patch_size+1;k<h1;k++)
    {
        for(int l=patch_size+1;l<w1;l++)
        {
            *(img+(l+k*h1)) = l+k*h1;
            //correlation_coefficient(k,l,diff)
        }
    }
     */

    //mexPrintf("Hello World!\n");
    return;
}