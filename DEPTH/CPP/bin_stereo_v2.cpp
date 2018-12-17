#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void map_epip_line(int map, int &c1, int &c2) {
    
    // horizontal epipolar line
    if (map == 1)
    {
        c1 = 0;
        c2 = 2;
    } 
    // 60 degrees (left top, right bottom)
    else if (map == 2)
    {
        c1 = 1;
        c2 = 1;
    }
    // 120 degrees (left bottom, right top)
    else if (map == 3)
    {
        c1 = 1;
        c2 = -1;
    }
}
    

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dimsL, *dimsR;
    int xL, yL, xR, yR;
    double * cost_volume, * imgL, * imgR;
    int i,j,x,y,h,w,d,dd;
  
    if (nrhs != 4)
    {
        mexErrMsgTxt("Four inputs required. Left, Center and Right Image (BW). Dmin, Dmax");
    }
    
    //figure out dimensions
    dimsL = mxGetDimensions(prhs[0]);
    dimsR = mxGetDimensions(prhs[1]);
    yL = (int)dimsL[0]; xL = (int)dimsL[1];
    yR = (int)dimsR[0]; xR = (int)dimsR[1];
    
    int dmin = (int)*mxGetPr(prhs[2]);
    int dmax = (int)*mxGetPr(prhs[3]);
    int numdisp = (dmax - dmin)+1;
    int pad = fmax(abs(dmin), abs(dmax));
    const mwSize cost_volume_dims[]={yL,xL,numdisp};
    plhs[0] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    imgL = mxGetPr(prhs[0]);
    imgR = mxGetPr(prhs[1]);
    cost_volume = mxGetPr(plhs[0]);
    mexPrintf("xL=%d, xR=%d\n", xL, xR);
    int ws = 5;
    int hws = floor(ws/2);
    
    double sad_sum = 0.0;
    double weight = 0.0;
    double diff = 0.0;
    double norm_sad = 0.0;
    
    for (x=pad+hws; x<xL-pad-hws; x++) // 199; x<200; x++) //
    {
        for (y=hws; y<yL-hws; y++) //199; y < 200; y++) //
        {
            for (d=dmin; d<=dmax; d++)
            {
                dd = d - dmin;
                //mexPrintf("x=%d, y=%d, dd=%d\n", x, y, dd);
                sad_sum = 0.0;
                for (h=-hws; h<=hws; h++)
                {
                    for (w=-hws; w<=hws; w++)
                    {
                        diff = fmin(50, std::abs(imgL[y+w +((x+h)*yL)] - imgR[y+w + ((x-d+h)*yL)]));
                        weight = exp(-sqrt(h*h+w*w) -diff/50) ;
                        sad_sum += weight * (diff*diff);
                    }
                }
                sad_sum /= (ws*ws);
                norm_sad = 1 - exp(-sad_sum/2.0);
                cost_volume[y + (x*yL) + (dd*xL*yL)] = norm_sad; //norm_sad;
                //cost_volume[y+w + ((x+h)*yL) + (dd*xL*yL)] = imgR[y+w +  ((x+h-d)*yL)];
            }
        }
    }
    
    
    return;
}
            