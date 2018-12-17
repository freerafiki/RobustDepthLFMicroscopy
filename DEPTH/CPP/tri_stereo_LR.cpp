#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

/**
 * It computes the cost volume cubes from three images in a row
 * It gives back the two costs for each matching 
 * LEFT-CENTER --> cost_volume_L_LR, cost_volume_L_RL
 * CENTER-RIGHT --> cost_volume_R_LR, cost_volume_R_RL
 * they can then be used for left-right consistency check and then fused
 */

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
    const mwSize *dimsL, *dimsR, *dimsC;
    int xL, yL, xR, yR, xC, yC;
    double * cost_volume_L, * cost_volume_R, * imgL, * imgR, *imgC;
    int i,j,k,h,w,d;
  
    if (nrhs != 8)
    {
        mexErrMsgTxt("Eight inputs required. Left, Center and Right Image (BW). Dmin, Dmax, Window Size, Alpha and Map");
    }
    
    //figure out dimensions
    dimsL = mxGetDimensions(prhs[0]);
    dimsC = mxGetDimensions(prhs[1]);
    dimsR = mxGetDimensions(prhs[2]);
    yL = (int)dimsL[0]; xL = (int)dimsL[1];
    yC = (int)dimsC[0]; xC = (int)dimsC[1];
    yR = (int)dimsR[0]; xR = (int)dimsR[1];
    
    int dmin = (int)*mxGetPr(prhs[3]);
    int dmax = (int)*mxGetPr(prhs[4]);
    int numdisp = (dmax - dmin)+1;
    const mwSize cost_volume_dims[]={yL,xL,numdisp};
    plhs[0] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    imgL = mxGetPr(prhs[0]);
    imgC = mxGetPr(prhs[1]);
    imgR = mxGetPr(prhs[2]);
    cost_volume_L_LR = mxGetPr(plhs[0]);
    cost_volume_L_RL = mxGetPr(plhs[1]);
    cost_volume_R_LR = mxGetPr(plhs[2]);
    cost_volume_R_RL = mxGetPr(plhs[3]);

    int ws = (int)*mxGetPr(prhs[5]);
    double alpha = *mxGetPr(prhs[6]);
    int hws = floor(ws/2);
    int dd = 0;
    int left_census[ws][ws];
    int central_census[ws][ws];
    int right_census[ws][ws];
    int map = *mxGetPr(prhs[7]);
    double sad_sumL = 0.0; double sad_sumR = 0.0;
    double TRUNC_AD = 10000; //50;
    double census_hdL = 0.0; double census_hdR = 0.0;
    double weightL = 0.0; double weightR = 0.0;
    double diffL = 0.0; double diffR = 0.0;
    const double THRESH_WEIGHT = 0.5;
    const double MAX_COST_SAD = ws*50.0;
    const double MAX_COST_CENSUS = ws*0.7;
    const double GAMMA_SAD = 50.0;
    const double GAMMA_CEN = 1.0;
    double norm_SAD_L = 0.0;
    double norm_CEN_L = 0.0;
    double norm_SAD_R = 0.0;
    double norm_CEN_R = 0.0;
    double g1 = (ws);
    double g2 = 10.0;
    int c1 = 0; int c2 = 0; int i_ = 0; int j_ = 0;
    double leftVal = 0.0; double rightVal = 0.0;
    int pad = std::max(std::abs(dmin), std::abs(dmax)); // we accept also negative disparities
    double weight_sumL = 0.0; double weight_sumR = 0.0;
    mexPrintf("dmin=%d, dmax=%d, ws=%d, hws=%d, alpha=%3.3f\n", dmin, dmax, ws, hws, alpha);
    for(i=pad+hws; i<xL-pad-hws; i++) //dmax+hws+1;i<s1-dmax-hws-1;i++) // 100;i<300;i++) // //(i=hws+1;i<s1-hws-1;i++)
    {
        for(j=pad+hws; j<yL-pad-hws; j++) //dmax+hws+1;j<s2-dmax-hws-1;j++) //  100;j<300;j++) //  (j=hws+1;j<s2-hws-1;j++)
        {       
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    if (imgC[j+h+(i+w)*yL] < imgC[j+(i)*yL])
                        central_census[h+hws][w+hws] = 0;
                    else
                        central_census[h+hws][w+hws] = 1;
                }
            } 
            //mexPrintf("\n\npixel %d, %d\n", i, j);
            // second round
            for (d=dmin; d<=dmax; d++) //dmax; d++)
            {   
                // d is the disparity, dd the index in the array (dd >= 0)
                dd = d - dmin; //numdisp - 1 - (d - dmin); // works for both positive and negative dmin;
                //mexPrintf("\nDisp=%d, index=%d ", d, dd);
                if (d == 0)
                {
                    sad_sumL = (MAX_COST_SAD); // /= (ws*TRUNC_AD);
                    sad_sumR = (MAX_COST_SAD);
                    census_hdL = (MAX_COST_CENSUS);
                    census_hdR = (MAX_COST_CENSUS);
                }
                else
                {
                    map_epip_line(map, c1, c2);
                    // takes into account epipolar line 
                    j_ = round(c1*0.866*(d));
                    i_ = round(c2*0.5*(d));
                    sad_sumL = 0.0;
                    census_hdL = 0.0;
                    sad_sumR = 0.0;
                    census_hdR = 0.0;
                    // using i_ and j_ for the two images (works for the 3 left and the 3 rights)
                    leftVal = imgL[j-j_+(i-i_)*yL];
                    rightVal = imgR[j+j_+(i+i_)*yL];
                    for(h=-hws;h<hws;h++)
                    {
                        for(w=-hws;w<hws;w++)
                        {
                            // build census windows
                            if (imgR[j+j_+h+(i+i_+w)*yL] < imgR[j+j_+(i+i_)*yL])
                                right_census[h+hws][w+hws] = 0;
                            else
                                right_census[h+hws][w+hws] = 1;
                            if (imgL[j+j_+h+(i-i_+w)*yL] < imgL[j+j_+(i-i_)*yL])
                                left_census[h+hws][w+hws] = 0;
                            else
                                left_census[h+hws][w+hws] = 1;

                            //weights --> using adaptive support windows
                            diffR = std::min(TRUNC_AD, std::abs(imgC[j+h+(i+w)*yL] - imgR[j+j_+h+(i+i_+w)*yL]));
                            diffL = std::min(TRUNC_AD, std::abs(imgC[j+h+(i+w)*yL] - imgL[j-j_+h+(i-i_+w)*yL]));
                            weightL = exp(-diffL / g2 -(sqrt((h*h)+(w*w))/g1));
                            weightR = exp(-sqrt((h*h)+(w*w))/g1 - diffR / g2);
                            //mexPrintf("weightL=%3.3f, weightR=%3.3f, h=%d, w=%d, LD=%3.3f, RD=%3.3f\n", weightL, weightR, h, w, std::abs(imgL[i+h+(j+w-d)*xL] - leftVal), std::abs(rightVal - imgR[i+h+(j+w+d)*xL]));
                            // add the weighted cost
                            census_hdL += weightL * std::abs(left_census[h+hws][w+hws] - central_census[h+hws][w+hws]);
                            census_hdR += weightR * std::abs(right_census[h+hws][w+hws] - central_census[h+hws][w+hws]);
                            sad_sumL += weightL * diffL;
                            sad_sumR += weightR * diffR;
                            weight_sumL += weightL;
                            weight_sumR += weightR;
                            //mexPrintf("actual diff=%3.3f", std::min(TRUNC_AD, std::abs(imgC[i+h+(j+w)*xL] - imgL[i-i_+h+(j-j_+w)*xL])));
                            //mexPrintf("LEFT: %3.3f - %3.3f = %3.3f\n", imgL[i+h+(j+w)*xL], imgL[i+h+(j+w-d)*xL], std::abs(imgL[i+h+(j+w)*xL] - imgL[i+h+(j+w-d)*xL]));
                            //mexPrintf("sad_sumL: %3.3f\n", sad_sumL);
                            //mexPrintf("censR=%d-%d=%d\n", right_census[h+hws][w+hws], central_census[h+hws][w+hws], std::abs(right_census[h+hws][w+hws] - central_census[h+hws][w+hws]));
                        }
                    }
                    //mexPrintf("before: sadL=%3.3f, sadR=%3.3f, censusL=%3.3f, censusR=%3.3f\n", sad_sumL, sad_sumR, census_hdL, census_hdR);
                    
                    // SHOULD THE SUM BE NORMALIZED ??
                    // sad_sumL /= (ws*ws); // divided by number of contributions
                    // sas_sumL /= weight_sumL; // divide by sum of the weights, so that weights sum up to 1!
                    sad_sumL = std::min(sad_sumL, MAX_COST_SAD); // /= (ws*TRUNC_AD);
                    sad_sumR = std::min(sad_sumR, MAX_COST_SAD);
                    census_hdL = std::min(census_hdL, MAX_COST_CENSUS);
                    census_hdR = std::min(census_hdR, MAX_COST_CENSUS);
                    norm_SAD_L = 1 - exp(-sad_sumL / GAMMA_SAD);
                    norm_SAD_R = 1 - exp(-sad_sumR / GAMMA_SAD);
                    norm_CEN_L = 1 - exp(-census_hdL / GAMMA_CEN);
                    norm_CEN_R = 1 - exp(-census_hdR / GAMMA_CEN);
                }

                //mexPrintf("after: sadL=%3.3f, sadR=%3.3f, censusL=%3.3f, censusR=%3.3f\n", sad_sumL, sad_sumR, census_hdL, census_hdR);
                //mexPrintf("normalized: sadL=%3.3f, sadR=%3.3f, censusL=%3.3f, censusR=%3.3f   ", norm_SAD_L, norm_SAD_R, norm_CEN_L, norm_CEN_R);
                //mexPrintf("cv=%3.3f\n", norm_SAD_L*alpha+norm_CEN_L*(1-alpha));
                cost_volume_L[j + i * yL + dd * (xL * yL)] = norm_SAD_L*alpha+norm_CEN_L*(1-alpha);
                cost_volume_R[j + i * yL + dd * (xL * yL)] = norm_SAD_R*alpha+norm_CEN_L*(1-alpha);
            }
        }
    }

    return;
}
            