# Bounding p-value for the null-hypothesis of (one-way) no signaling in a Bell test

Matlab codes to compute a p-value upper bound for the null-hypothesis of one-way or two-way no-signaling based on the raw data obtained in a Bell test:
- [Min_KLDivergence_NSCond.m](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/Min_KLDivergence_NSCond.m): For a given bipartite relative frequency F={f(a,b|x,y)} and input distribution Pxy(x,y), find the (unconditional) correlation minimizing the Kullback–Leibler divergence from F to the one-way or two-way no-signaling sets
- [MaxBellValue_NSCond.m](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/MaxBellValue_NSCond.m): For a given bipartite Bell expression specified by {beta(a,b,x,y)}, find the maximal Bell value allowed by correlations satisfying either the one-way or two-way no-signaling conditions
- [PBRPValueBoundForNSFromRawData.m](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/PBRPValueBoundForNSFromRawData.m): Take the data obtained in a bipartite Bell test and compute a p-value upper bound from the adapted prediction-based-ratio (PBR) protocol.

As an example, the Matlab data file [Data_Washington_91_98_CircNL_BellTest15.mat](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/Data_Washington_91_98_CircNL_BellTest15.mat) contains the complete set of data for one of the Bell tests performed on IBMQ Washington for qubits 91 and 98.  

The codes require
- [YALMIP](https://yalmip.github.io/) - a MATLAB toolbox for optimization modeling
- [MOSEK](https://www.mosek.com/) - a commercial convex optimization solver

These codes have been tested on Matlab R2020b, Yalmip version 22-June-2023, and Mosek version 9.3.
