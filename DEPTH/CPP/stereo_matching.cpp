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
        mexErrMsgTxt("Two inputs required. Left and Right Image (BW), dmin and dmax");
    }
    
    //figure out dimensions
    dimsL = mxGetDimensions(prhs[0]);
    dimsR = mxGetDimensions(prhs[1]);
    yL = (int)dimsL[0]; xL = (int)dimsL[1];
    yR = (int)dimsR[0]; xR = (int)dimsR[1];
    
    int dmin = *mxGetPr(prhs[2]);
    int dmax = *mxGetPr(prhs[3]);
    int numdisp = dmax - dmin;
    const mwSize cost_volume_dims[]={xL,yL,numdisp};
    plhs[0] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(3, cost_volume_dims, mxDOUBLE_CLASS, mxREAL);
    imgL = mxGetPr(prhs[0]);
    imgR = mxGetPr(prhs[1]);
    cost_volume_sad = mxGetPr(plhs[0]);
    cost_volume_census = mxGetPr(plhs[1]);
    cost_volume_ncc = mxGetPr(plhs[2]);
    cost_volume_ssd = mxGetPr(plhs[3]);
    int ws = 9;
    int hws = floor(ws/2);
    int dd = 0;
    int left_census[ws][ws];
    int right_census[ws][ws];
    double sad_sum = 0.0;
    double TRUNC_AD = 50;
    double ssd_sum = 0.0;
    double TRUNC_SD = 1000;
    double census_hd = 0.0;
    double mean_p1 = 0.0; double mean_p2 = 0.0; double mean_p12 = 0.0; double sigma_p1 = 0.0;
    double sigma_p2 = 0.0; double tmp_p1 = 0.0; double tmp_p2 = 0.0;
    double weight = 0.0;
    double g1 = 1.0;
    double g2 = 3.0;
    double leftVal = 0.0; double rightVal = 0.0;
    int pad = std::max(std::abs(dmin), std::abs(dmax)); // we accept also negative disparities
    double weight_sum = 0.0;
    mexPrintf("x=%d, y=%d, dmin=%d, dmax=%d", xL, yL, dmin, dmax);
    for(i=180;i<200;i++) //pad+hws; i<xL-pad-hws; i++) //dmax+hws+1;i<s1-dmax-hws-1;i++) //  //(i=hws+1;i<s1-hws-1;i++)
    {
        for(j=300;j<310;j++) // pad+hws; j<yL-pad-hws; j++) //dmax+hws+1;j<s2-dmax-hws-1;j++) //   (j=hws+1;j<s2-hws-1;j++)
        {
            // mean of the window
            mean_p1 = 0.0;
            mean_p2 = 0.0;            
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    if (imgL[i+h+(j+w)*xL] < imgL[i+(j)*xL])
                        left_census[h+hws][w+hws] = 0;
                    else
                        left_census[h+hws][w+hws] = 1;
                    mean_p1 += imgL[(i+h)*xL+(j+w)];
                    mean_p2 += imgR[(i+h)*xL+(j+w)];
                }
            } 
            mean_p1 /= (ws*ws);
            mean_p2 /= (ws*ws);
            
            // second round

            for (d=dmin; d<dmax; d++)
            {
                sigma_p1 = 0.0;
                sigma_p2 = 0.0;
                mean_p12 = 0.0;
                sad_sum = 0.0;
                ssd_sum = 0.0;
                census_hd = 0.0;
                for(h=-hws;h<hws;h++)
                {
                    for(w=-hws;w<hws;w++)
                    {
                        if (imgR[i+h+(j+d+w)*xL] < imgR[i+(j+d)*xL])
                            right_census[h+hws][w+hws] = 0;
                        else
                            right_census[h+hws][w+hws] = 1;
                        leftVal = imgL[i+h+(j+w)*xL];
                        rightVal = imgR[i+h+(j+w+d)*xL];
                        //if (std::abs(leftVal - rightVal) < 0.5)
                        //    weight = 1.0;
                        //else
                        weight = exp(-sqrt((h*h)+(w*w))/g1 + std::abs(imgL[i+h+(j+w)*xL] - imgL[i+(j)*xL]) / g2);
                        //exp(-sqrt((h*h)+(w*w))) ; 
                        // ;
                        census_hd += weight * std::abs(right_census[h+hws][w+hws] - left_census[h+hws][w+hws]);
                        ssd_sum += weight * std::min(TRUNC_SD, (imgL[i+h+(j+w)*xL] - imgR[i+h+(j+w+d)*xL])*(imgL[i+h+(j+w)*xL] - imgR[i+h+(j+w+d)*xL]));
                        sad_sum += weight * std::min(TRUNC_AD, std::abs(imgL[i+h+(j+w)*xL] - imgR[i+h+(j+w+d)*xL]));
                        weight_sum += weight;
                        mexPrintf("c=%1.5f, s=%1.5f, a=%1.5f, w=%1.5f, h=%d, w=%d\n", weight * std::abs(right_census[h+hws][w+hws] - left_census[h+hws][w+hws]), weight * std::min(TRUNC_SD, (leftVal - rightVal)*(imgL[i+h+(j+w)*xL] - imgR[i+h+(j+w+d)*xL])), std::min(TRUNC_AD, std::abs(imgL[i+h+(j+w)*xL] - imgR[i+h+(j+w+d)*xL])), weight, h, w);
                        //tmp_p1 = (imgL[(i+h)*xL+(j+w)] - mean_p1);
                        //tmp_p2 = (imgR[(i+h)*xL+(j+w)] - mean_p2);
                        //mean_p12 += tmp_p1*tmp_p2;
                        //sigma_p1 += pow(tmp_p1,2);
                        //sigma_p2 += pow(tmp_p2,2);
                    }
                }
                //sigma_p1 /= (ws*ws);
                //sigma_p2 /= (ws*ws);
                //mean_p12 /= (ws*ws);
                sad_sum /= (ws*ws*TRUNC_AD*weight_sum);
                ssd_sum /= (ws*ws*TRUNC_SD*weight_sum);
                census_hd /= (ws*ws*weight_sum);
                //mexPrintf("sad=%3.3f, ssd=%3.3f, census=%3.3f\n", sad_sum ,ssd_sum, census_hd);
                // d is the disparity, dd the index in the array (dd >= 0)
                dd = d - dmin; // works for both positive and negative dmin;
                cost_volume_sad[i + j * xL + dd * (xL * yL)] = sad_sum;
                cost_volume_census[i + j * xL + dd * (xL * yL)] = census_hd;
                //cost_volume_ncc[i + j * xL + dd * (xL * yL)] = (mean_p12) / (sigma_p1 * sigma_p2);
                cost_volume_ssd[i + j * xL + dd * (xL * yL)] = ssd_sum;
            }
        }
    }

    return;
}
            