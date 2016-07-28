#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    int node, f, s,s1,s2,feature, num_states,num_nodes,num_cols,num_nodesW,num_features,
            *x,*y,*featureStart; /* alternatively mwsize */
    
    double O,*unit,*gr_start,*gr_end,*bin;
    
    unit = mxGetPr(prhs[0]);
    gr_start = mxGetPr(prhs[1]);
    gr_end = mxGetPr(prhs[2]);
    bin = mxGetPr(prhs[3]);
    x = (int*)mxGetPr(prhs[4]);
    y = (int*)mxGetPr(prhs[5]);
    num_states = ((int*)mxGetPr(prhs[6]))[0];
    featureStart = (int*)mxGetPr(prhs[7]);
    
    if (!mxIsClass(prhs[5],"int32"))
        mexErrMsgTxt("y must be int32");
    
    num_nodesW = mxGetDimensions(prhs[0])[0];
    num_nodes = mxGetDimensions(prhs[4])[0];
    num_features = mxGetDimensions(prhs[4])[1];
    
    for(node = 0;node < num_nodes;node++) {
        for(f = 0;f < num_features;f++) {
            if (x[node + num_nodes*f] != 0) {
                /*printf("n = %d, f = %d\n",n,f);*/
                for(s = 0;s < num_states;s++) {
                    if (s == y[node]-1)
                        O = 1;
                    else
                        O = 0;                         
                    unit[featureStart[f]-1 + x[node + num_nodes*f]-1 + num_nodesW*s] += O;                   
                }
            }
        }
    }
    
    for(s = 0; s < num_states; s++) {
        if (s == y[0]-1)
            O = 1;
        else
            O = 0;
        gr_start[s] += O;
        if (s == y[num_nodes-1]-1)
            O = 1;
        else
            O = 0;
        gr_end[s] += O;
    }
    
    for(node = 0; node < num_nodes-1; node++) {
        for(s1 = 0; s1 < num_states; s1++) {
            for(s2 = 0; s2 < num_states; s2++) {
                if (s1 == y[node]-1 && s2 == y[node+1]-1)
                    O = 1;
                else
                    O = 0;
                bin[s1 + num_states*s2] += O;
            }
        }
    }
    
}
