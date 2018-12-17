#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dims, *dims_fl;
    int x, y, fl1, fl2;
    double * output, *input, *values;
    int i,j,h,w;
  
    if (nrhs != 2)
    {
        mexErrMsgTxt("Two inputs required.");
    }
    
    // depth map
    dims = mxGetDimensions(prhs[0]);
    x = (int)dims[0]; y = (int)dims[1];
    dims_fl = mxGetDimensions(prhs[1]);
    fl1 = (int)dims_fl[0];
    
    //mexPrintf("h1=[%d], w1=[%d], h2=[%d], w2=[%d]\n", dimx1, numdims1, dimx2, numdims2);
    
    plhs[0] = mxCreateDoubleMatrix(x, y, mxREAL);
    output = mxGetPr(plhs[0]);
    input = mxGetPr(prhs[0]);
    values = mxGetPr(prhs[1]);
    //do something
    for(i=0;i<x;i++)
    {
        for(j=0;j<y;j++)
        {
            output[i*y+j] = values[(int)(input[i*y+j])];
        }
    }
    return;
}