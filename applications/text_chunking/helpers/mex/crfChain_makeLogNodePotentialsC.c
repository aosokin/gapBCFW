#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n,f,state,featureParam,
            *x,*featureStart,
            num_states,num_nodes,num_features,num_featuresTotal,num_rows;
    double *wv,*logNodePot;
    
    x = (int*)mxGetPr(prhs[0]);
    wv = mxGetPr(prhs[1]);
    featureStart = (int*)mxGetPr(prhs[2]);
    num_states = ((int*)mxGetPr(prhs[3]))[0];
    
    if (!mxIsClass(prhs[0],"int32")||!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32"))
        mexErrMsgTxt("{x,featureStart,num_states} must be int32");
    
    num_rows = mxGetDimensions(prhs[0])[0];
    num_nodes = num_rows;
    num_features = mxGetDimensions(prhs[0])[1];
    num_featuresTotal = featureStart[num_features]-1;
    
    plhs[0] = mxCreateDoubleMatrix(num_nodes,num_states,mxREAL);
    logNodePot = mxGetData(plhs[0]);
    
    for(n = 0; n < num_nodes; n++) {
        for(f = 0; f < num_features; f++) {
            if (x[n + num_rows*f] != 0) {
                for(state = 0; state < num_states; state++) {
                    featureParam = featureStart[f]-1 + x[n + num_rows*f]-1;
                    logNodePot[n + num_nodes*state] += wv[featureParam + num_featuresTotal*state];
                }
            }
        }
    }
    
    for(state = 0; state < num_states; state++) {
        logNodePot[num_nodes*state] += wv[num_featuresTotal*num_states + state];
        logNodePot[num_nodes-1 + num_nodes*state] += wv[num_featuresTotal*num_states + num_states + state];
    }
}
