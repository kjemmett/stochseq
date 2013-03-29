#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get pointers to inputs
    double *px_z, *p;
    px_z = mxGetPr(prhs[0]);
    p = mxGetPr(prhs[1]);
    
    //double *px_z, *A, *At, *p;
    //px_z = mxGetPr(prhs[0]);
    //A = mxGetPr(prhs[1]);
    //At = mxGetPr(prhs[2]);
    //p = mxGetPr(prhs[3]);

    //mwIndex *A_ir, *At_ir, *A_jc, *At_jc;
    //A_ir = mxGetIr(prhs[1]);
    //A_jc = mxGetJc(prhs[1]);
    //At_ir = mxGetIr(prhs[2]);
    //At_jc = mxGetJc(prhs[2]);

    // loop indices
    long i, j, k, l, t;

	// dimensions: T, K
    // K=L
	mwSize T, K;
	T = mxGetM(prhs[0]);
    K = mxGetN(prhs[0]);

	// // local variables: a, b, c
	float *a = new float[T*K];
	float *b = new float[T*K];
	float *c = new float[T]; 

    // local variables: pr, pf
    // transition probabilities
    double pf = *p;
    //mexPrintf("hello %.2f\n", pf);
    double pr = 1 - *p;

	// initialize to zero
	for (i=0; i<T*K; i++)
	{
		a[i] = 0;
		b[i] = 0;
	}
	//for (i=0; i<T; i++)
	//{
	//	c[i] = 0;
	//}

	// Forward Sweep - Calculate
	//
	//   a(t, k)  =  sum_l px_z(t,k) A(l, k) alpha(t-1, l)  
	//   c(t)     =  sum_k a(t, k)
	//
	// and normalize 
	//
	//   a(t, k)  /=  c(t)

	//// a(0, k)  =  px_z(0, k) a1(k)
	//// normalize a(0,k) by c(k)

    // we can explicitly initialize alpha(T=0)
    // and directly set a[0], c[0]
    a[0] = 1;
    c[0] = px_z[0];

    for (t = 1; t < T; t++)
    {
        // a(t, k)  =  sum_l px_z(t,k) A(l, k) alpha(t-1, l)  
        
        // transition into first state
        a[t] = px_z[t] * pr * a[T+t-1];
        c[t] = a[t];
        // transition into second state
        a[T+t] = px_z[T+t] * (a[t-1] + pr*a[2*T + t-1]);
        c[t] += a[T+t];
        // transition into states 3 to K-1
        for (k = 2; k < K-1; k++) {
            a[k*T+t] = px_z[k*T+t] * (pr*a[(k+1)*T+t-1] + pf*a[(k-1)*T+t-1]);
            c[t] += a[k*T + t];
        }
        // transition into final state K
        a[(K-1)*T+t] = px_z[(K-1)*T+t] * pf*a[(k-1)*T+t-1];
        c[t] += a[(K-1)*T+t];

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

    b[K*T-1] = 1;
    d = a[K*T-1] * b[K*T-1];
    b[K*T-1] /= d;
    
	// t = T-2:0
	for (t = T-2; t >= 0; t--)
	{
		// b(t, k)  =  sum_l px_z(t+1,l) A(k, l) beta(t+1, l)  
        b[t] = px_z[T+t+1] * b[T+t+1] / c[t+1];
        for (k=1; k < K-1; k++) {
            b[k*T+t] = (px_z[(k+1)*T+t+1]*pf*b[(k+1)*T+t+1] + px_z[(k-1)*T+t+1]*pr*b[(k-1)*T+t+1]) / c[t+1];
        }
        b[(K-1)*T+t] = px_z[(K-2)*T+t+1]*pr*b[(K-2)*T+t+1] / c[t+1];
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
