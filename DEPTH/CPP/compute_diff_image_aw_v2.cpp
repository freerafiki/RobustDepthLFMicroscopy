#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dimsAIF, *dimsCFS;
    int dimxA, dimyA, dimxF, dimyF;
    double * output, *aif, *cur_fs;
    int i,j,h,w;
  
    if (nrhs != 4)
    {
        mexErrMsgTxt("Four inputs required. AIF Image, FS Image, Window Size and Alpha.");
    }

    //inputs
    dimsAIF = mxGetDimensions(prhs[0]);
    dimsCFS = mxGetDimensions(prhs[1]);
    dimyA = (int)dimsAIF[0]; dimxA = (int)dimsAIF[1];
    dimyF = (int)dimsCFS[0]; dimxF = (int)dimsCFS[1];
    aif = mxGetPr(prhs[0]);
    cur_fs = mxGetPr(prhs[1]);
    int ws = (int)*mxGetPr(prhs[2]);
    double alpha = *mxGetPr(prhs[3]);
    int hws = floor(ws/2.0);
    //outputs
    plhs[0] = mxCreateDoubleMatrix(dimyA, dimxA, mxREAL);
    output = mxGetPr(plhs[0]);
    mexPrintf("y=%d, x=%d", dimyA, dimxA);
    
    // COMPUTATIONS
    // NCC and SAD variables
    double mean_p1 = 0.0;
    double mean_p2 = 0.0;
    double mean_p12 = 0.0;
    double sigma_p1 = 0.0;
    double sigma_p2 = 0.0;
    double ncc_val = 0.0;
    double tmp_p1 = 0.0;
    double tmp_p2 = 0.0;
    double sad_sum = 0.0;
    double sad_val = 0.0;
    double ncc_val_low = 0.0;
    
    const double THRESH_WEIGHT = 1.0;
    const double MAX_COST_SAD = 100.0;
    const double MIN_COST_NCC = 0.0;
    const double GAMMA_SAD = 20.0;
    const double GAMMA_NCC = 0.5;
    double norm_SAD = 0.0;
    double norm_NCC = 0.0;
    
    double TRUNC_AD = 50; // truncate AD
    double g1 = (ws*ws); // sigma for distances
    double g2 = 25.0; // sigma for intensity (if range changes (not 0-255 but 0-1) change accordingly)
                        // --> (image[0,255] -> g2 = 25.0 +- 10 ; image[0,1] -> g2 = 0.25 +- 0.1)
    double weight = 0.0;
    double weight_sum = 0.0; // for normalization purposes
    
    //mexPrintf("ws=%d, alpha=%3.3f", ws, alpha);
    //LOOP
    for(i=hws+1;i<dimxA-hws-1;i++) // 350; i<400; i++) //
    {
        for(j=hws+1;j<dimyA-hws-1;j++) // 350; j<400; j++) //
        {
            // mean of the window
            mean_p1 = 0.0;
            mean_p2 = 0.0;
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {
                    mean_p1 += aif[(j+w)*dimyA+(i+h)];
                    mean_p2 += cur_fs[(j+w)*dimyA+(i+h)];
                }
            }
            mean_p1 /= (ws*ws);
            mean_p2 /= (ws*ws);
            
            // real calculations
            sigma_p1 = 0.0;
            sigma_p2 = 0.0;
            mean_p12 = 0.0;
            sad_sum = 0.0;
            for(h=-hws;h<hws;h++)
            {
                for(w=-hws;w<hws;w++)
                {

                    weight = exp(-(sqrt((h*h)+(w*w))/g1 + std::abs(aif[(j+w)*dimyA+(i+h)] - cur_fs[(j+w)*dimyA+(i+h)]) / g2));
                    
                    tmp_p1 = (aif[(i+h)*dimyA+(j+w)] - mean_p1);
                    tmp_p2 = (cur_fs[(i+h)*dimyA+(j+w)] - mean_p2);
                    mean_p12 += weight * tmp_p1*tmp_p2;
                    sigma_p1 += weight * pow(tmp_p1,2);
                    sigma_p2 += weight * pow(tmp_p2,2);
                    sad_sum += weight * std::min(TRUNC_AD, std::abs(aif[(j+w)*dimyA+(i+h)] - cur_fs[(j+w)*dimyA+(i+h)]));
                    weight_sum += weight;
                }
            }   
            sigma_p1 /= (ws*ws);
            sigma_p2 /= (ws*ws);
            mean_p12 /= (ws*ws);
            sad_sum /= (ws*ws);
            weight_sum /= (ws*ws);
            
            //ncc has large values for correct matches
            ncc_val = std::max(MIN_COST_NCC, (mean_p12) / sqrt(sigma_p1 * sigma_p2));
            ncc_val_low = 1 - ncc_val;
            //mexPrintf("ncc_val = %3.3f\n, meanp12 = %3.3f, sigma_p1 = %3.3f, sigma_p2 = %3.3f",ncc_val, mean_p12, sigma_p1, sigma_p2);
            //sad has low values for correct matches
            sad_val = std::min(MAX_COST_SAD, sad_sum);
            
            norm_SAD = 1 - exp(-sad_val / GAMMA_SAD);
            norm_NCC = 1 - exp(-ncc_val_low / GAMMA_NCC);
            //mexPrintf("normSAD=%3.3f, normNCC=%3.3f", norm_SAD, norm_NCC);
            
            //mexPrintf("sigma_cur[%f], sigma_aif[%f], sigma_cross[%f] \n",sigma_cur, sigma_aif, sigma_cross);
            output[(j)*dimyA+(i)] = (alpha)*norm_NCC + (1-alpha)*norm_SAD; //-cur_fs[i*dimy1+j]; 
        }
    }
    
}