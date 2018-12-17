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
    double sum = 0.0;
    //do something
    for(i=hws+1;i<dimx1-hws-1;i++)
    {
        for(j=hws+1;j<dimy1-hws-1;j++)
        {
            sum = 0.0;

            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {

                    sum += (aif[(i+h)*dimy1+(j+w)] - cur_fs[(i+h)*dimy1+(j+w)])*(aif[(i+h)*dimy1+(j+w)] - cur_fs[(i+h)*dimy1+(j+w)]);
                }
            }   
            //mexPrintf("sigma_cur[%f], sigma_aif[%f], sigma_cross[%f] \n",sigma_cur, sigma_aif, sigma_cross);
            output[i*dimy1+j] = sum / (ws*ws); //-cur_fs[i*dimy1+j]; 
        }
    }

    return;
}