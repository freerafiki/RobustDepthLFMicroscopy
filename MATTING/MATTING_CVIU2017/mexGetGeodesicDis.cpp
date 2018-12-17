#include "mex.h"
#include "math.h"

double findMinD( double *I, double *D, double *J, int m, int n, int ak[8], int i, int j)
{
	double minD, gradient, gamma, temp;
	int k,x,y, N;

	N = m*n;
	minD = D[ i + j*m ];
	for( k = 1; k <= 4; k++)
	{
		x = j + ak[2*k-2];
		y = i + ak[2*k-1];
		if( y < 0 || x < 0 || y >= m || x >= n )
			continue;
		else
		{
			gradient =  pow( I[ i + j*m ] - I[ y + x*m], 2) + pow( I[ i+j*m+N ] - I[ y+x*m+N ], 2) + pow( I[ i+j*m+2*N ] -I[ y+x*m+2*N], 2);
			gamma = 1/(J[i + j*m] + 0.01);
			temp = D[y + x*m] + pow( ak[2*k-2]*ak[2*k-2] + ak[2*k-1]*ak[2*k-1] + gamma*gamma*gradient, 0.5);
			if( temp < minD)
				minD = temp;
		}
	}
	return minD;
}

void rasterScan ( double *I, double *D, double *J, int m, int n)
{
	int Ak[4][8] = {{-1,-1,-1,0,-1,1,0,-1},{1,1,1,0,1,-1,0,1},{1,-1,0,-1,-1,-1,1,0},{-1,1,0,1,1,1,-1,0}};
	int i,j,k,ak[8];
	double minD;

	for(k = 0; k <= 7; k++)
		ak[k] = Ak[0][k];
	for(j = 0; j < n; j++)
		for(i = 0; i < m; i++)
		{
			D[ i + j*m ] = findMinD(I,D,J,m,n,ak,i,j);
		};

	for(k = 0; k <= 7; k++)
		ak[k] = Ak[1][k];
	for(j = n-1; j >= 0; j--)
		for(i = m-1; i >= 0; i--)
		{
			D[ i + j*m ] = findMinD(I,D,J,m,n,ak,i,j);
		};

	for(k = 0; k <= 7; k++)
		ak[k] = Ak[2][k];
	for(i = 0; i <= m-1; i++)
		for(j = n-1; j >= 0; j--)
		{
			D[ i + j*m ] = findMinD(I,D,J,m,n,ak,i,j);
		};

	for(k = 0; k <= 7; k++)
		ak[k] = Ak[3][k];
	for(i = m-1; i >= 0; i--)
		for(j = 0; j <= n-1; j++)
		{
			D[ i + j*m ] = findMinD(I,D,J,m,n,ak,i,j);
		};

}




void mexFunction ( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	double *I, *J, *D, *M;
	int m,n,gamma,k,v;

	I = mxGetPr(prhs[0]);
	J = mxGetPr(prhs[1]);
	M = mxGetPr(prhs[2]);
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	v = 10000;

	
	plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);

	D = mxGetPr(plhs[0]);
	for( k = 0; k <= m*n -1; k++)
		D[k] = v*M[k];

	rasterScan(I,D,J,m,n);

}
