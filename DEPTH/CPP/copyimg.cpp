#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dimsL, *dimsR;
    int xL, yL, xR, yR;
    double * output, *img1, *img2;
    int i,j,h,w,d;
    if (nrhs != 2)
    {
        mexErrMsgTxt("Two inputs required.");
    }
    
    //first is cost_volume
    dimsL = mxGetDimensions(prhs[0]);
    dimsR = mxGetDimensions(prhs[1]);
    yL = (int)dimsL[0]; xL = (int)dimsL[1];
    yR = (int)dimsR[0]; xR = (int)dimsR[1];
    int ws = 15;
    int hws = floor(ws/2);
    double gs = 0.0;
    plhs[0] = mxCreateDoubleMatrix(yL, xL, mxREAL);
    img1 = mxGetPr(prhs[0]);
    img2 = mxGetPr(prhs[1]);
    output = mxGetPr(plhs[0]);
    mexPrintf("xL=%d, yL=%d\nxR=%d, yR=%d\n",xL, yL,xR,yR);
    for(i=ws;i<xL-ws;i++) //(j=hws+1;j<s2-hws-1;j++)   
    {
        for(j=ws;j<yL-ws;j++) //(i=hws+1;i<s1-hws-1;i++)
        {
            gs = 0.0;
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    gs += img1[j+w + (i)*yL];
                }
            }
            gs /= (ws*ws);
            output[j + i * yL] = gs; // 1 / sum;
        }
    }
    

    return;
}