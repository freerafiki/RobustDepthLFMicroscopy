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
    double * cost_volume_sad, * cost_volume_census, * cost_volume_ncc, * cost_volume_ssd, * imgL, * imgR;
    int i,j,k,h,w,d;
  
    if (nrhs != 4)
    {
        mexErrMsgTxt("Four inputs required. Left and Right Image (BW), dmin and dmax");
    }
    
    //figure out dimensions
    dimsL = mxGetDimensions(prhs[0]);
    dimsR = mxGetDimensions(prhs[1]);
    yL = (int)dimsL[0]; xL = (int)dimsL[1];
    yR = (int)dimsR[0]; xR = (int)dimsR[1];
    
    int dmin = *mxGetPr(prhs[2]);
    int dmax = *mxGetPr(prhs[3]);
    int numdisp = dmax - dmin;
    const mwSize cost_volume_dims[]={yL,xL,numdisp};
    plhs[0] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    imgL = mxGetPr(prhs[0]);
    imgR = mxGetPr(prhs[1]);
    cost_volume_sad = mxGetPr(plhs[0]);
    cost_volume_census = mxGetPr(plhs[1]);
    int ws = 21;
    int hws = floor(ws/2);
    int dd = 0;
    double sad_sum = 0.0;
    double census_hd = 0.0;
    double weight = 0.0;
    int left_census[ws][ws];
    int right_census[ws][ws];
    int pad = std::max(std::abs(dmin), std::abs(dmax));
    mexPrintf("x=%d, y=%d, dmin=%d, dmax=%d", xL, yL, dmin, dmax);
    for(i=pad+hws+10; i<xL-pad-hws-10; i++) //dmax+hws+1;i<s1-dmax-hws-1;i++) // 180;i<200;i++) // //(i=hws+1;i<s1-hws-1;i++)
    {
        for(j=pad+hws+10; j<yL-pad-hws-10; j++) //dmax+hws+1;j<s2-dmax-hws-1;j++) // 300;j<310;j++) //  (j=hws+1;j<s2-hws-1;j++)
        {
            //mexPrintf("i=%d, j=%d\n", i,j);
            for (d=dmin; d<dmax; d++)
            {
                for(h=-hws;h<hws;h++)
                {
                    for(w=-hws;w<hws;w++)
                    {
                        if (imgL[j+h+(i+w)*yL] < imgL[j+(i)*yL])
                            left_census[h+hws][w+hws] = 0;
                        else
                            left_census[h+hws][w+hws] = 1;
                    }
                } 
                sad_sum = 0.0;
                census_hd = 0.0;
                for(h=-hws;h<hws;h++)
                {
                    for(w=-hws;w<hws;w++)
                    {
                        if (imgR[j+h+(i+w+d)*yR] < imgR[j+(i+d)*yR])
                            right_census[h+hws][w+hws] = 0;
                        else
                            right_census[h+hws][w+hws] = 1;
                        weight = 1; //exp(-sqrt(h*h+w*w)/10 - std::abs(imgL[j+w+(i+h)*yL] - imgR[j+w+(i+h+d)*yL])/50);
                        census_hd += weight * std::abs(right_census[h+hws][w+hws] - left_census[h+hws][w+hws]);
                        sad_sum += weight * (std::abs(imgL[j+w+(i+h)*yL] - imgR[j+w+(i+h+d)*yL]));
                        //mexPrintf("c=%1.5f, s=%1.5f, a=%1.5f, w=%1.5f, h=%d, w=%d\n", weight * std::abs(right_census[h+hws][w+hws] - left_census[h+hws][w+hws]), weight * std::min(TRUNC_SD, (leftVal - rightVal)*(imgL[i+h+(j+w)*xL] - imgR[i+h+(j+w+d)*xL])), std::min(TRUNC_AD, std::abs(imgL[i+h+(j+w)*xL] - imgR[i+h+(j+w+d)*xL])), weight, h, w);
                    }
                }
                sad_sum /= (ws*ws);
                //sad_sum /= (ws*ws); //*weight_sum);
                //ssd_sum /= (ws*ws*weight_sum);
                census_hd /= (ws*ws); //*weight_sum);
                //mexPrintf("sad=%3.3f, ssd=%3.3f, census=%3.3f\n", sad_sum ,ssd_sum, census_hd);
                // d is the disparity, dd the index in the array (dd >= 0)
                dd = d - dmin; // works for both positive and negative dmin;
                cost_volume_sad[j + i * yL + dd * (xL * yL)] = sad_sum;
                cost_volume_census[j + i * yL + dd * (xL * yL)] = census_hd;
            }
        }
    }

    return;
}
            