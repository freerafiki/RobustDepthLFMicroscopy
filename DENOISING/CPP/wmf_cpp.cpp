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
    plhs[0] = mxCreateDoubleMatrix(dimxI, dimyI, mxREAL);
    output = mxGetPr(plhs[0]);
    
    // COMPUTATIONS
    double g1 = ws; // sigma for distances
    double g2 = 250.0; // sigma for intensity (if range changes (not 0-255 but 0-1) change accordingly)
                        // --> (image[0,255] -> g2 = 25.0 +- 10 ; image[0,1] -> g2 = 0.25 +- 0.1)
    double weight_median = 0.0;
    double weight_sum = 0.0; // for normalization purposes
    double weight_limit = 0.0;
    double weight = 0.0;
    int vec_size = ws*ws;
    int row[vec_size];
    int col[vec_size];
    double weights[vec_size];
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
            for(h=-hws;h<=hws;h++)
            {
                for(w=-hws;w<=hws;w++)
                {
                    row_ind = h+hws;
                    col_ind = w+hws;
                    row[(row_ind)*ws+col_ind] = i+h;
                    col[(row_ind)*ws+col_ind] = j+w;
                    weight =  exp(-(sqrt((h*h)+(w*w))/g1)); // + std::abs(ref[(i+h)*dimyI+(j+w)] - ref[(i)*dimyI+(j)]) / g2));
                    //mexPrintf("r=%d, c=%d, i=%d - ",row_ind, col_ind, row_ind*ws+col_ind);
                    if (weight < 1 & weight > 0)
                    {
                        weights[row_ind*ws+col_ind] = weight;
                        weight_sum += weight;
                    }
                    else
                    {
                        weights[row_ind*ws+col_ind] = 0.0;
                        weight_sum += 0.0;
                    }
                }
            }   
            // select the median value
            //mexPrintf("medianweight = %3.3f, weight_sum=%3.3f, weigth_limit=%3.3f\n", weight_median, weight_sum, weight_limit);
            weight_median = weight_sum / 2.0; // central value
            weight_sum = 0.0;
            weight_limit = 0.0; // here we sum every time to see the cell we are in
            v = 0;
            while (!found)
            {
                weight_limit += weights[v]; // adding up each value
                mexPrintf("w-%dth = %3.3f", v, weight_limit);
                // when we pass the central value we arrived to our central value
                if (v >= (vec_size))
                {
                    chosen_x = -1;
                    chosen_y = -1;
                    found = true;
                }
                if (weight_limit > weight_median)
                {
                    // arrived
                    chosen_x = row[v];
                    chosen_y = col[v];
                    found = true;
                    mexPrintf("x=%d,y=%d,chosen_x=%d,chosen_y=%d\n", i,j, chosen_x, chosen_y);
                }
                else
                {   
                    v++; 
                }
            }
            
            if (chosen_x > -1)
            {
                med_val = img[chosen_x*dimyI+chosen_y];
                mexPrintf("value = %3.3f, old val = %3.3f, x=%d, y=%d, \n", med_val, img[i*dimyI+j], i, j);
            }
            else
            {   
                med_val = img[i*dimyI+j];
                mexPrintf("w_lim = %3.3f, x=%d, y=%d, value = %3.3f\n", weight_limit, i, j, med_val);
            }
            found = false;
            //mexPrintf("medval = %3.3f, chosen_x=%d, chosen_y=%d\n", med_val, chosen_x, chosen_y);
            output[i*dimyI+j] = med_val; //-cur_fs[i*dimy1+j]; 
        }
    }
    
}