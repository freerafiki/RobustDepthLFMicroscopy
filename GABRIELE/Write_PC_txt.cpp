#define _USE_MATH_DEFINES
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <cmath>
#include<stdlib.h>
#include<stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //variables
    const mwSize *dimsPC_map, *dimyEI_name;
    int dimyX, dimyY, dimyZ, dimyRGBZ;
    int h,w,v,cont,r,g,b,dimyEI_r,dimyEI_c;
    double x,y,z,correzX,correzY;
    double *PC_map, *EI_naming;
    
    FILE *fd;
    
    //inputs
    dimsPC_map = mxGetDimensions(prhs[0]);
    dimyX = (int)dimsPC_map[0]; dimyY = (int)dimsPC_map[1];
    dimyZ = (int)dimsPC_map[2]; dimyRGBZ = (int)dimsPC_map[3];
    
    PC_map = mxGetPr(prhs[0]);
    EI_naming = mxGetPr(prhs[1]);
    dimyEI_name = mxGetDimensions(prhs[1]);
    dimyEI_r = (int)dimyEI_name[0]; dimyEI_c = (int)dimyEI_name[1]; 
    mexPrintf("r0 = %d", dimyEI_r);
    cont=0;
    fd = fopen("PC_4.txt", "w");

    for (int h = 0; h < dimyEI_r; h++)
    {
        for (int w = 0; w < dimyY; w++)
        {
            for (int v = 0; v < dimyX; v++)
            {   
                correzY = round((2-EI_naming[h])*cos(M_PI/6)*(PC_map[w+v*dimyX+4*h*dimyX*dimyY+3*dimyX*dimyY]));
                correzX = round((2-EI_naming[dimyEI_r+h])*0.5*(PC_map[w+v*dimyX+4*h*dimyX*dimyY+3*dimyX*dimyY]));
                y = w+correzY;
                x = v+correzX;
                z = 10*PC_map[w+v*dimyX+4*h*dimyX*dimyY+3*dimyX*dimyY];
                r = (int)PC_map[w+v*dimyX+4*h*dimyX*dimyY+0*dimyX*dimyY];
                g = (int)PC_map[w+v*dimyX+4*h*dimyX*dimyY+1*dimyX*dimyY];
                b = (int)PC_map[w+v*dimyX+4*h*dimyX*dimyY+2*dimyX*dimyY];
                
                fprintf(fd, "%.2f\t%.2f\t%.2f\t%d\t%d\t%d\n",
                        (float)x,(float)y,(float)z,//(float)x,(float)y,(float)z,
                        (int)r,
                        (int)g,
                        (int)b);

            }
        }
    }
    
fclose(fd);
mexPrintf("%d",cont);
return;
}