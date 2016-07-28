#define INF 1E20
#include <math.h>
#include <sys/types.h>
#include "mex.h"
#include <iostream>
/*
 * distance_transform.cpp
 * this code is adapted from [1]. However we adapted it to be able to deal with 
 * undealt situations where the parabolas were linear functions (a_c and a_p equal 0).
 *
 * [1] Articulated Pose Estimation by a Graphical Model with Image Dependent
 * Pairwise Relations, Chen and Yuille
 *
 * distance transforms algorithm for quadratic function ax^2 + bx, based on Felzenswalb and Huttenlocher.
 * @article{felzenszwalb2005pictorial,
 *   title={Pictorial structures for object recognition},
 *   author={Felzenszwalb, Pedro F and Huttenlocher, Daniel P},
 *   journal={International Journal of Computer Vision (IJCV)},
 *   year={2005},
 * }
 */

static inline double square(double x) { return x*x; }

void dt1d(double *src, double *dst, int *ptr, int step, int len,
        double a_c, double b_c, double a_p, double b_p, double dshift_c, double dshift_p, int dlen) {
  int   *v = new int[len];
  float *z = new float[len+1];
  int k = 0;
  int q = 0;
  v[0] = 0;
  z[0] = -INF;
  z[1] = +INF;
  
 // special case where parabolas become linear function
  if ((a_c>-1e-7) && (a_p>-1e-7)){
	double max_value = -INF;
	int val_k = -1;
	 
	for (k = 0; k <= dlen-1; k++) {
		double s = k*(b_c-b_p)-b_c*dshift_c-b_p*dshift_p+src[k*step];
		if (s>max_value) {
			val_k = k;
		    max_value = s;
		}
	}

	for (q = 0; q <= dlen-1; q++) {
		dst[q*step] = src[val_k*step] + a_c * square(q + dshift_c - val_k) + b_c * -(q + dshift_c - val_k)
		+ a_p * square(q - dshift_p - val_k) + b_p * (q - dshift_p - val_k);
		ptr[q*step] = val_k;
	  }

  } else{
	  for (q = 1; q <= len-1; q++) {
		float s = ( (src[q*step] - src[v[k]*step])
		- b_c * -(q - v[k]) + a_c * (square(q) - square(v[k]))
		- b_p * (q - v[k]) + a_p * (square(q) - square(v[k]))
		+ 2*a_c * (q-v[k])*(-dshift_c) + 2*a_p * (q-v[k])*(dshift_p) )
		/ ( 2*a_c*(q-v[k]) + 2*a_p*(q-v[k]));
		
		while (s <= z[k]) {
		  k--;
		  s = ( (src[q*step] - src[v[k]*step])
		  - b_c * -(q - v[k]) + a_c * (square(q) - square(v[k]))
		  - b_p * (q - v[k]) + a_p * (square(q) - square(v[k]))
		  + 2*a_c * (q-v[k])*(-dshift_c) + 2*a_p * (q-v[k])*(dshift_p) )
		  / ( 2*a_c*(q-v[k]) + 2*a_p*(q-v[k]));
		}
		k++;
		v[k]   = q;
		z[k]   = s;
		z[k+1] = +INF;
	  }
	  
	  k = 0;
	  for (q = 0; q <= dlen-1; q++) {
		while (z[k+1] < q)
		  k++;
		dst[q*step] = src[v[k]*step] + a_c * square(q + dshift_c - v[k]) + b_c * -(q + dshift_c - v[k])
		+ a_p * square(q - dshift_p - v[k]) + b_p * (q - dshift_p - v[k]);
		ptr[q*step] = v[k];
	  }
  }
  delete [] v;
  delete [] z;
}



// matlab entry point
// [M, Ix, Iy] = distance_transform(vals, ax, bx, ay, by)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 9)
    mexErrMsgTxt("Wrong number of inputs");
  if (nlhs != 3)
    mexErrMsgTxt("Wrong number of outputs");
  
  // ---- values
  double *vals = (double *)mxGetPr(prhs[0]);
  if (mxGetNumberOfDimensions(prhs[0]) != 2
          || mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Invalid input val");
  
  int sizx  = mxGetN(prhs[0]);
  int sizy  = mxGetM(prhs[0]);
  
  // ---- deformation weight child->parent
  if (mxGetNumberOfElements(prhs[1]) != 4
          || mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("Invalid child->parent deformation weights");
  }
  double *dw_c = mxGetPr(prhs[1]);
  double ax_c = -dw_c[0];            // 2nd order
  double bx_c = -dw_c[1];            // 1st order
  double ay_c = -dw_c[2];
  double by_c = -dw_c[3];
  
  // ---- deformation weight parent->child
  if (mxGetNumberOfElements(prhs[2]) != 4
          || mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("Invalid parent->child deformation weights");
  }
  double *dw_p = mxGetPr(prhs[2]);
  double ax_p = -dw_p[0];
  double bx_p = -dw_p[1];
  double ay_p = -dw_p[2];
  double by_p = -dw_p[3];
  
  // ---- child->parent mean
  if (mxGetNumberOfElements(prhs[3]) != 2
          || mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("Invalid child->parent mean");
  }
  double *mean_c = mxGetPr(prhs[3]);
  
  // ---- child->parent var
  if (mxGetNumberOfElements(prhs[4]) != 2
          || mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("Invalid child->parent var");
  }
  double *var_c = mxGetPr(prhs[4]);
  
  // ---- parent->child mean
  if (mxGetNumberOfElements(prhs[5]) != 2
          || mxGetClassID(prhs[5]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("Invalid parent->child mean");
  }
  double *mean_p = mxGetPr(prhs[5]);
  
  // ---- parent->child var
  if (mxGetNumberOfElements(prhs[6]) != 2
          || mxGetClassID(prhs[6]) != mxDOUBLE_CLASS) {
    mexErrMsgTxt("Invalid parent->child var");
  }
  double *var_p = mxGetPr(prhs[6]);
  
  int32_t lenx  = (int32_t)mxGetScalar(prhs[7]);
  int32_t leny  = (int32_t)mxGetScalar(prhs[8]);
  
  mxArray  *mxM = mxCreateNumericMatrix(leny,lenx,mxDOUBLE_CLASS, mxREAL);
  mxArray *mxIy = mxCreateNumericMatrix(leny,lenx,mxINT32_CLASS, mxREAL);
  mxArray *mxIx = mxCreateNumericMatrix(leny,lenx,mxINT32_CLASS, mxREAL);
  double   *M = (double *)mxGetPr(mxM);
  int32_t *Iy = (int32_t *)mxGetPr(mxIy);
  int32_t *Ix = (int32_t *)mxGetPr(mxIx);
  
  double   *tmpM =  (double *)mxCalloc(leny*sizx, sizeof(double));
  int32_t *tmpIy = (int32_t *)mxCalloc(leny*sizx, sizeof(int32_t));
  
  for (int x = 0; x < sizx; x++)
    dt1d(vals+x*sizy, tmpM+x*leny, tmpIy+x*leny, 1, sizy,
            ay_c/square(var_c[1]), by_c/var_c[1],
            ay_p/square(var_p[1]), by_p/var_p[1],
            mean_c[1], mean_p[1], leny);
  
  for (int y = 0; y < leny; y++)
    dt1d(tmpM+y, M+y, Ix+y, leny, sizx,
            ax_c/square(var_c[0]), bx_c/var_c[0],
            ax_p/square(var_p[0]), bx_p/var_p[0],
            mean_c[0], mean_p[0], lenx);
  
  // get argmins and adjust for matlab indexing from 1
  for (int x = 0; x < lenx; x++) {
    for (int y = 0; y < leny; y++) {
      int p = x*leny+y;
      Iy[p] = tmpIy[Ix[p]*leny+y] + 1;
      Ix[p] = Ix[p] + 1;
    }
  }
  
  mxFree(tmpM);
  mxFree(tmpIy);
  plhs[0] = mxM;
  plhs[1] = mxIx;
  plhs[2] = mxIy;
  return;
}
