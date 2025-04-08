# Bounding p-value for the null-hypothesis of (one-way) no signaling in a Bell test

Matlab codes to compute a p-value upper bound for the null-hypothesis of one-way or two-way no-signaling based on the raw data obtained in a Bell test:
- [Min_KLDivergence_NSCond.m](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/Min_KLDivergence_NSCond.m): For a given bipartite relative frequency F={f(a,b|x,y)} and input distribution Pxy(x,y), find the (unconditional) correlation minimizing the Kullback–Leibler divergence from F to the one-way or two-way no-signaling sets
- [MaxBellValue_NSCond.m](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/MaxBellValue_NSCond.m): For a given bipartite Bell expression specified by {beta(a,b,x,y)}, find the maximal Bell value allowed by correlations satisfying either the one-way or two-way no-signaling conditions
- [PBRPValueBoundForNSFromRawData.m](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/PBRPValueBoundForNSFromRawData.m): Take the data obtained in a bipartite Bell test and compute a p-value upper bound from the adapted prediction-based-ratio (PBR) protocol.

As an example, the Matlab data file [Data_Washington_91_98_CircNL_BellTest15.mat](https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/Data_Washington_91_98_CircNL_BellTest15.mat) contains the complete set of data for one of the Bell tests performed on IBMQ Washington for qubits 91 and 98. After loading the Matlab data file, one can execute the code PBRPValueBoundForNSFromRawData.m by feeding the loaded "BellTestData" to obtain the p-value bound.

The complete data set used to obtain the results presented in https://arxiv.org/abs/2503.18949 is given in the zip file [IBMQ_Data.zip]( https://github.com/ycliangTW/AlmostDICertCrossTalks/blob/main/IBMQ_Data.zip). Upon unzipping the file, one finds five folders labeled by the name of the IBM device. Inside each folder, for each qubit pair tested, there are a pair of files containing, respectively, the input data (x,y) and the output data (a,b). The raw data are listed as an integer from 0 to 3 where 0,1,2,3 are the decimals for the corresponding two bits: 0 = 00, 1 = 01, 2 = 10, 3 = 11. Each column in the output data file consists of all the output data for one complete Bell test consisting of 1800 trials.

The codes require
- [YALMIP](https://yalmip.github.io/) - a MATLAB toolbox for optimization modeling
- [MOSEK](https://www.mosek.com/) - a commercial convex optimization solver

These codes have been tested on Matlab R2020b, Yalmip version 22-June-2023, and Mosek version 9.3.
