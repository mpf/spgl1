#include "mex.h"

/* ----------------------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* ----------------------------------------------------------------------- */
{  mxArray       *vectorX;
   mxArray       *vectorY;
   double        *ptrX, *ptrY, *ptrSqrt1, *ptrSqrt2, x, t;
   int            i, d, transpose;

   /* Free memory and exit if no parameters are given */
   if (nrhs == 0)
   {  if (nlhs != 0) { mexErrMsgTxt("No output arguments expected."); }
      return ;
   }

   /* Check for proper number of arguments */
   if (nrhs != 4) { mexErrMsgTxt("Four input arguments required."); }
   if (nlhs >  1) { mexErrMsgTxt("Too many output arguments.");             }

    /* Extract the arguments */
    vectorX   = (mxArray *)prhs[0];
    transpose = (mxGetScalar(prhs[1]) == 0 ? 0 : 1);
    ptrSqrt1  = mxGetPr((mxArray *)prhs[2]);
    ptrSqrt2  = mxGetPr((mxArray *)prhs[3]);

    if (!mxIsDouble(vectorX)  || ((mxGetM(vectorX) > 1) &&
       (mxGetN(vectorX) > 1)) || (mxGetNumberOfDimensions(vectorX) != 2))
    {   mexErrMsgTxt("Parameter 'x' has to be a double vector.");
    }

    /* Get the problem size */
    d = mxGetNumberOfElements(vectorX);

    /* Create output array */
    if (transpose == 1)
    {  vectorY = mxCreateDoubleMatrix(d-1, 1, mxREAL);
       ptrX = mxGetPr(vectorX);
       ptrY = mxGetPr(vectorY);
       
       t = 0; x = *ptrX;
       ptrY --; ptrSqrt1 --; ptrSqrt2 --; /* We index Y and sqrt starting from 1 */
       for (i = 0; i < d-1; )
       {  i ++;
          t -= x;
          x  = ptrX[i];
          ptrY[i] = ptrSqrt1[i] * t + ptrSqrt2[i] * x;
       }
    }
    else
    {  vectorY = mxCreateDoubleMatrix(d+1, 1, mxREAL);
       ptrX = mxGetPr(vectorX);
       ptrY = mxGetPr(vectorY);
       
       /* ptrSqrt1: sqrt(1 / (i * (i+1))) */
       /* ptrSqrt2: sqrt(i / (i + 1.0)) */
       
       t = 0;
       ptrX --; ptrSqrt1 --; ptrSqrt2 --; /* We index X starting from 1 */
       for (i = d; i > 0; i--)
       {  x = ptrX[i];
          ptrY[i] = t + ptrSqrt2[i] * x;
          t -= ptrSqrt1[i] * x;
       }
       ptrY[0] = t;
    }

    /* Set result */
    plhs[0] = vectorY;
    
    return ;
}
