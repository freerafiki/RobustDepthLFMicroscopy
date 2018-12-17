#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dimsIMG, *dimsREF;
    int dimxI, dimyI, dimxR, dimyR;
    double *output, *img, *ref;
    int i,j,h,w,v;
  
    if (nrhs != 3)
    {
        mexErrMsgTxt("Three inputs required. Image, Reference Image, Window Size.");
    }

    //inputs
    dimsIMG = mxGetDimensions(prhs[0]);
    dimsREF = mxGetDimensions(prhs[1]);
    dimyI = (int)dimsIMG[0]; dimxI = (int)dimsIMG[1];
    dimyR = (int)dimsREF[0]; dimxR = (int)dimsREF[1];
    img = mxGetPr(prhs[0]);
    ref = mxGetPr(prhs[1]);
    int ws = (int)*mxGetPr(prhs[2]);
    int hws = floor(ws/2.0);
    //outputs
    plhs[0] = mxCreateDoubleMatrix(dimyI, dimxI, mxREAL);
    output = mxGetPr(plhs[0]);
    
    // COMPUTATIONS
    double g1 = 10; // sigma for distances
    double g2 = 25; // sigma for intensity (if range changes (not 0-255 but 0-1) change accordingly)
                        // --> (image[0,255] -> g2 = 25.0 +- 10 ; image[0,1] -> g2 = 0.25 +- 0.1)
    double bf_sum = 0.0;
    double weight_sum = 0.0; // for normalization purposes
    double weight = 0.0;
    int vec_size = ws*ws;
    int row_ind = 0;
    int col_ind = 0;
    int chosen_x = 0;
    int chosen_y = 0;
    bool found = false;
    double med_val = 0.0;
    //mexPrintf("ws=%d, alpha=%3.3f", ws, alpha);
    //LOOP
    for(i= hws+1;i<dimxI-hws-1;i++)  //50; i<450; i++) //
    {
        for(j= hws+1;j<dimyI-hws-1;j++)  //50; j<450; j++) //
        {
            //mexPrintf("\n");
            bf_sum = 0;
            weight_sum = 0;
            for(h=-hws;h<=hws;h++)
            {
                for(w=-hws;w<=hws;w++)
                {
                    weight = 1-exp(-(sqrt((h*h)+(w*w))/g1)-(ref[(i+h)*dimyI+(j+w)] - ref[(i)*dimyI+(j)])*(ref[(i+h)*dimyI+(j+w)] - ref[(i)*dimyI+(j)]) / g2); // -(sqrt((h*h)+(w*w))/g1) 
                    bf_sum += img[(i+h)*dimyI+(j+w)]*weight;
                    weight_sum += weight;
                }
            }   
            
            bf_sum /= (weight_sum);
            //mexPrintf("medval = %3.3f, chosen_x=%d, chosen_y=%d\n", med_val, chosen_x, chosen_y);
            output[i*dimyI+j] = bf_sum; //-cur_fs[i*dimy1+j]; 
        }
    }
    
}