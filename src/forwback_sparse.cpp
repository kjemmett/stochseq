#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get pointers to inputs
    double *px_z, *A, *At, *a1, *bT;
    px_z = mxGetPr(prhs[0]);
    A = mxGetPr(prhs[1]);
    At = mxGetPr(prhs[2]);
    a1 =  mxGetPr(prhs[3]);
    bT =  mxGetPr(prhs[4]);

    mwIndex *A_ir, *At_ir, *A_jc, *At_jc;
    A_ir = mxGetIr(prhs[1]);
    A_jc = mxGetJc(prhs[1]);
    At_ir = mxGetIr(prhs[2]);
    At_jc = mxGetJc(prhs[2]);

    // loop indices
    long i, j, k, l, t;

	// dimensions: T, K
	mwSize T, K;
	T = mxGetM(prhs[0]);
    K = mxGetN(prhs[0]);

	// // local variables: a, b, c
	float *a = new float[T*K];
	float *b = new float[T*K];
	float *c = new float[T]; 

    // create outputs for a, b, c
    // plhs[2] = mxCreateDoubleMatrix(T, K, mxREAL);
    // plhs[3] = mxCreateDoubleMatrix(T, K, mxREAL);
    // plhs[4] = mxCreateDoubleMatrix(T, 1, mxREAL);
    // double *a, *b, *c;
    // a = mxGetPr(plhs[2]);
    // b = mxGetPr(plhs[3]);
    // c = mxGetPr(plhs[4]);

	// initialize to zero
	for (i=0; i<T*K; i++)
	{
		a[i] = 0;
		b[i] = 0;
	}
	for (i=0; i<T; i++)
	{
		c[i] = 0;
	}

	// Forward Sweep - Calculate
	//
	//   a(t, k)  =  sum_l px_z(t,k) A(l, k) alpha(t-1, l)  
	//   c(t)     =  sum_k a(t, k)
	//
	// and normalize 
	//
	//   a(t, k)  /=  c(t)

	// a(0, k)  =  px_z(0, k) a1(k)
	for (k = 0; k < K; k++) 
	{
		//a[k*T] = a1[k] * px_z[k*T];
	    a[k*T] = a1[k];
    	c[0] += a[k*T];
	}
	// normalize a(0,k) by c(k)
	for (k = 0; k < K; k++) 
	{
		a[k*T] /= c[0];
	}

    for (t = 1; t < T; t++)
    {
        // a(t, k)  =  sum_l px_z(t,k) A(l, k) alpha(t-1, l)  
        for (k = 0; k < K; k++) 
        {
            // sparse loop over row indices for column in A(:, k)
            for (i = A_jc[k]; i < A_jc[k+1]; i++)  
            {
                // A(l, k) == A[i]
                l = A_ir[i];
                // a(t,k) += px_z(t,k) A(l, k) alpha(t-1, l)  
                a[k*T + t] += px_z[k*T + t] * A[i] * a[l*T + t-1];
            }
            // c(t) += a(t,k)
            c[t] += a[k*T + t];
        }
        // normalize a(t,k) by c(t)
        for (k = 0; k < K; k++) 
        {
            a[k*T + t] /= c[t];
        }
    }
		
	// Back sweep - calculate
	//
	// b(t,k)  =  1/c(t+1) sum_l px_z(t+1, l) A(k, l) beta(t+1, l) 

	// b(T-1,k) = 1
    float d = 0;
	for (k = 0; k < K; k++) 
	{
		b[k*T + (T-1)] = bT[k];
        d += a[k*T + (T-1)] * b[k*T + (T-1)];
	}
    for (k = 0; k < K; k++) {
        b[k*T + (T-1)] /= d;
    }

	// t = T-2:0
	for (t = T-2; t >= 0; t--)
	{
		// b(t, k)  =  sum_l px_z(t+1,l) A(k, l) beta(t+1, l)  
        for (k = 0; k < K; k++) 
        {
            // sparse loop over column indices A(k, :)
            for (j = At_jc[k]; j < At_jc[k+1]; j++)  
            {
                // A(k, l) == At[j]
                l = At_ir[j];
                // b(t ,k) += px_z(t+1, l) A(k, l) betal(t+1, l)  
                b[k*T + t] += px_z[l*T + t+1] * At[j] * b[l*T + t+1];
            }          
            // normalize b(t,k) by c(t+1)
            b[k*T + t] /= c[t+1];
        }
	}

	// allocate outputs: g, lnZ
	// plhs[0] = mxCreateNumericMatrix(T, K, mxSINGLE_CLASS, mxREAL);
	// float *g = (float *) mxGetData(plhs[0]);
    plhs[0] = mxCreateDoubleMatrix(T, K, mxREAL);
    double *g = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *lnZ = mxGetPr(plhs[1]);

	// g(t,k) = a(t,k) * b(t,k)
	for (i=0; i<T*K; i++)
	{
		g[i] = a[i] * b[i];
	}

	// ln_Z = sum_t log(c[t])
	lnZ[0] = 0;
	for (t=0; t<T; t++)
	{
		lnZ[0] += log(c[t]);
	}

    free(a);
    free(b);
    free(c);

	return;
}
