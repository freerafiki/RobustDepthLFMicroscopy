#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dims1, *dims2;
    int dimx1, dimy1, numdims1;
    double * output, *input;
    int i,j,h,w;
  
    if (nrhs != 3)
    {
        mexErrMsgTxt("Three inputs required.");
    }
    
    //figure out dimensions
    dims1 = mxGetDimensions(prhs[0]);
    numdims1 = mxGetNumberOfDimensions(prhs[0]);
    dimy1 = (int)dims1[0]; dimx1 = (int)dims1[1];
    int step_pix = *mxGetPr(prhs[1]);
    int offset = *mxGetPr(prhs[2]);
    
    //mexPrintf("h1=[%d], w1=[%d], h2=[%d], w2=[%d]\n", dimx1, numdims1, dimx2, numdims2);
    
    plhs[0] = mxCreateDoubleMatrix(dimx1, dimy1, mxREAL);
    output = mxGetPr(plhs[0]);
    input = mxGetPr(prhs[0]);
    //do something
    for(i=0;i<dimx1;i++)
    {
        for(j=0;j<dimy1;j++)
        {
            output[i*dimy1+j] = offset - step_pix*input[i*dimy1+j];
        }
    }
    return;
}