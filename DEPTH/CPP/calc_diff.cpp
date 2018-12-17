#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dims1, *dims2, *dims3, *dims4, *dims5;
    int dimx1, dimy1, dimx2, dimy2, dimx3, dimy3, dimx4, dimy4, dimx5, dimy5;
    double * output_sp1, *output_sp2, *output_points, *aif, *cur_fs, *sp1, *sp2, *mask;
    int i,j,h,w, x,y, k;
  
    if (nrhs != 7)
    {
        mexErrMsgTxt("Seven inputs required. \n1) Focal stack image \n2) all-in-focus image \n3) superpixels1 \n4) superpixels2 \n5) mask \n6) window size \n7) alpha\n");
    }
    
    // !!! BOTH IMAGES NEED TO BE GRAYSCALE!!! 
    // first input is focal stack image
    cur_fs = mxGetPr(prhs[0]);
    dims1 = mxGetDimensions(prhs[0]);
    dimy1 = (int)dims1[0]; dimx1 = (int)dims1[1];
    // second input is all-in-focus image
    aif = mxGetPr(prhs[1]);
    dims2 = mxGetDimensions(prhs[1]);
    dimy2 = (int)dims2[0]; dimx2 = (int)dims2[1];
    // third are the big superpixels
    sp1 = mxGetPr(prhs[2]);
    dims3 = mxGetDimensions(prhs[2]);
    dimy3 = (int)dims3[0]; dimx3 = (int)dims3[1];
    // fourth are the big superpixels
    sp2 = mxGetPr(prhs[3]);
    dims4 = mxGetDimensions(prhs[3]);
    dimy4 = (int)dims4[0]; dimx4 = (int)dims4[1];
    // fifth is the mask
    mask = mxGetPr(prhs[4]);
    dims5 = mxGetDimensions(prhs[4]);
    dimy5 = (int)dims5[0]; dimx5 = (int)dims5[1];
    // sixth is window size
    int ws = *mxGetPr(prhs[5]);
    int hws = floor(ws/2.0);
    // seventh alpha
    double alpha = *mxGetPr(prhs[6]);;
    double beta = 1 - alpha;
    // maximum superpixel 1
    int max_sp1 = sp1[dimy1*(dimy1-1)+dimx1-1];
    // maximum superpixel 2
    int max_sp2 = sp2[dimy1*(dimy1-1)+dimx1-1];
    
    // output will be three matrices with the costs
    plhs[0] = mxCreateDoubleMatrix(dimx1, dimy1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(dimx1, dimy1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(dimx1, dimy1, mxREAL);
    output_sp1 = mxGetPr(plhs[0]);
    output_sp2 = mxGetPr(plhs[1]);
    output_points = mxGetPr(plhs[2]);

    // variables to calculate ncc
    double ncc_sum = 0.0;
    double mean_aif = 0.0;
    double mean_cur = 0.0;
    double sigma_aif = 0.0;
    double sigma_cur = 0.0;
    double sigma_cross = 0.0;
    double cur_cur = 0.0;
    double cur_aif = 0.0;
    double sad_sum = 0.0;
    
    // to fill the superpixels parts
    double sp1_values [max_sp1];
    double sp2_values [max_sp2];
    int sp1_counter [max_sp1];
    int sp2_counter [max_sp2];
    int sp1_valid [max_sp1];
    int sp2_valid [max_sp2];
    
    //loop on the image
    for(i=hws+1;i<dimx1-hws-1;i++) // 300; i < 350; i++) //
    {
        for(j=hws+1;j<dimy1-hws-1;j++) //300; j < 350; j++) // 
        {
            if (mask[i*dimy1+j] > 0)
            {
                mean_aif = 0.0;
                mean_cur = 0.0;
                for(h=-hws;h<hws;h++)
                {
                    for(w=-hws;w<hws;w++)
                    {
                        mean_aif += aif[(i+h)*dimy1+(j+w)];
                        mean_cur += cur_fs[(i+h)*dimy1+(j+w)];
                    }
                }
                mean_aif /= (ws*ws);
                mean_cur /= (ws*ws);
                ncc_sum = 0.0;
                sigma_cur = 0.0;
                sigma_aif = 0.0;
                sigma_cross = 0.0;
                sad_sum = 0.0;
                //mexPrintf("mean_aif[%f], mean_cur[%f] \n",mean_aif, mean_cur);
                //mexPrintf("sigma_cur[%f], sigma_aif[%f], sigma_cross[%f] \n",sigma_cur, sigma_aif, sigma_cross);
                for(h=-hws;h<hws;h++)
                {
                    for(w=-hws;w<hws;w++)
                    {
                        cur_cur = (pow(cur_fs[(i+h)*dimy1+(j+w)] - mean_cur,2));
                        sigma_cur += cur_cur;
                        cur_aif = (pow(aif[(i+h)*dimy1+(j+w)] - mean_aif,2));
                        sigma_aif += cur_aif;
                        sigma_cross += (cur_cur)*(cur_aif);
                        sad_sum += std::abs(cur_fs[(i+h)*dimy1+(j+w)] - aif[(i+h)*dimy1+(j+w)]);
                    }
                }   
                ncc_sum = (sigma_aif*sigma_cur) / sigma_cross;
                //mexPrintf("sigma_cur[%f], sigma_aif[%f], sigma_cross[%f], sum=%.3f \n",sigma_cur, sigma_aif, sigma_cross, sum);
                if (sigma_cross == 0)
                    ncc_sum = 1.0;
                output_points[i*dimy1+j] = alpha*ncc_sum + (1-alpha)*sad_sum; //-cur_fs[i*dimy1+j];
                // fill the superpixels
                // first round is to fill the array with the values
                sp1_values[(int)(sp1[i*dimy1+j])] += (alpha*ncc_sum + (1-alpha)*sad_sum);
                sp2_values[(int)(sp2[i*dimy1+j])] += (alpha*ncc_sum + (1-alpha)*sad_sum);
                sp1_valid[(int)(sp1[i*dimy1+j])] += 1 ;
                sp2_valid[(int)(sp2[i*dimy1+j])] += 1 ;
            }
            sp1_counter[(int)(sp1[i*dimy1+j])] += 1 ;
            sp2_counter[(int)(sp2[i*dimy1+j])] += 1 ;
        }
    }
    // second round take from array and put into image
    //mexPrintf("second round");
    //for (k=0;k<max_sp1;k++)
    //    mexPrintf("sp1_values[%d]=%.3f\n",k,sp1_values[k]);
    int index_sp1 = 0;
    int index_sp2 = 0;
    for(i=hws+1;i<dimx1-hws-1;i++) // 100; i < 150; i++) //
    {
        for(j=hws+1;j<dimy1-hws-1;j++) // 300; j < 350; j++) //
        {
            index_sp1 = (int)(sp1[i*dimy1+j]);
            index_sp2 = (int)(sp2[i*dimy1+j]);
            if (sp1_valid[index_sp1] > (sp1_counter[index_sp1] / 2))
                output_sp1[i*dimy1+j] = sp1_values[index_sp1] / sp1_counter[index_sp1];
            if (sp2_valid[index_sp2] > (sp2_counter[index_sp2] / 2))
                output_sp2[i*dimy1+j] = sp2_values[index_sp2] / sp2_counter[index_sp2];
        }
    }
    
    
    return;
}