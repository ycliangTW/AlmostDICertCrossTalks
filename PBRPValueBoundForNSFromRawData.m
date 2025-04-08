% PBRPValueBoundForNSFromRawData
% requires: yalmip, mosek, MaxBellValue_NSCond, Min_KLDivergence_NSCond

% author: Yeong-Cherng Liang
%
%    [pValNS, RNS, DiagMinKL, DiagMaxBell] = PBRPValueBoundForNSFromRawData(BellTestData, Dir, Pxy, BellScenario, NumEst, pCutoff)
%    computes the p-value upper bound for the one-way or two-way no-
%    signaling (NS) hypothesis based on the adapted prediction-based-ratio
%    (PBR) protocol and "BellTestData" collected in a Bell test defined 
%    in the "BellScenario" and inputs generated according to the 
%    distribution "Pxy". The actual null hypothesis to be tested is 
%    specified by "Dir", taking values
%      'NS' for two-way NS, 'A2B' for one-way NS from A to B, and 'B2A' for
%       one-way NS from B to A
%    The number of trials used to evaluate the PBR is specified by "NumEst". 
%    "BellTestData" should be a table consistsing of 4 columns such that
%    the i-th row consists of the numbers [a_i b_i x_i y_i], i.e, the
%    outputs and inputs of the i-th experimental trial. "BellScenario"
%    should be an array of the form [nA nB nX nY] indicating the 
%    cardinality of the number of outpus for each input of both parties 
%    followed their number of inputs (we only consider a symmetrical Bell 
%    scenario where each input consists of the same number of outpus. The 
%    input distribution Pxy(x,y) is again a 2-dimensional array. "pCutoff"
%    specifies the p-value cutoff for the rejection of the null hypothesis.
%
%    The output of the code "pValNS" is the computed p-value upper bound,
%    "RNS" is the PBR calculated from the relative frequencies (from the 
%    first "NumEst" estimation trials) and the KL divergence computation.
%    "DiagMinKL" contains diagnostics information of the KL divergence
%    minimization, and "DiagMaxBell" contains diagnostics information
%    arising from Bell value maximization over the set of NS correlations 
%    specified by "Dir".

%% Copyright (C) 2025 Yeong-Cherng Liang, last modified on 7 Apr 2025

function [pValNS, RNS, DiagMinKL, DiagMaxBell] = PBRPValueBoundForNSFromRawData(BellTestData, Dir, Pxy, BellScenario, NumEst, pCutoff)



% Validate input
if ~ismatrix(BellTestData) 
    error('BellTestData must be a 2D array.');
else
    [NTrials Sz] = size(BellTestData);
    if Sz~=4 || length(BellScenario)~=4
        error('The code only works for a bipartite Bell test');
    end
end

% Extract the number of outputs
nA = BellScenario(1);
nB = BellScenario(2);

% Extract the number of inputs
nX = BellScenario(3);
nY = BellScenario(4);

options = sdpsettings;
options.verbose=0;
options.cachesolvers=1;

PWhite = ones(nA,nB,nX,nY)/4;

% Counts and frequencies for all data in a Bell test
NXY = zeros(nX,nY);
NABXY = zeros(nA,nB,nX,nY);
Freq = zeros(nA,nB,nX,nY);

% Counts and frequencies for the estimation/ training data
NXYEst = zeros(nX,nY);
NABXYEst = zeros(nA,nB,nX,nY);
FreqEst = zeros(nA,nB,nX,nY);

% Counts and frequencies for the test data
NXYTest = zeros(nX,nY);
NABXYTest = zeros(nA,nB,nX,nY);
FreqTest = zeros(nA,nB,nX,nY);

for x=1:nX
    for y=1:nY
        

        Idx = find(BellTestData(:,3)==x & BellTestData(:,4)==y);
        IdxEst = find(BellTestData(1:NumEst,3)==x & BellTestData(1:NumEst,4)==y);
        NXY(x,y) = length(Idx);
        NXYEst(x,y) = length(IdxEst);
        
        NXYTest(x,y) = NXY(x,y)-NXYEst(x,y);
        
        for a=1:nA
            for b=1:nB
                
                NABXY(a,b,x,y) = length(find(BellTestData(Idx,1)==a & BellTestData(Idx,2)==b));
                NABXYEst(a,b,x,y) = length(find(BellTestData(IdxEst,1)==a & BellTestData(IdxEst,2)==b));
                NABXYTest(a,b,x,y) = NABXY(a,b,x,y)-NABXYEst(a,b,x,y);
                
                Freq(a,b,x,y) = NABXY(a,b,x,y)/NXY(x,y);
                FreqEst(a,b,x,y) = NABXYEst(a,b,x,y)/NXYEst(x,y);
                FreqTest(a,b,x,y) = NABXYTest(a,b,x,y)/NXYTest(x,y);
                
            end
        end
        
    end
end

% If the relative frequency contains zero entries, we mix it with a bit of 
% white noise s.t. the numerics don't run into any unnecessary complication

if min(FreqEst(:))==0
    FreqEst = FreqEst*NumEst/(NumEst+1) + PWhite/(NumEst+1);
end

% Minimization of the KL divergence from the target set to the relative 
% frequencies: output "Q_opt" is the optimized unconditional distribution
if strcmp(Dir,'NS')
    [Q_opt, DiagMinKL] = Min_KLDivergence_NSCond(FreqEst, Pxy);
else
    [Q_opt, DiagMinKL] = Min_KLDivergence_NSCond(FreqEst, Pxy,Dir);
end

if DiagMinKL.problem~=0
    fprintf('Problem encountered during the KL divergence optimization!\n');
end

% Initialization of the PBR and the Bell expression
RNS = zeros(nA,nB,nX,nY);
betaNS = zeros(nA,nB,nX,nY);

for x=1:nX
    for y=1:nY
        for a=1:nA
            for b=1:nB
                
                % Ratio of the relative frequency to the optimizing
                % conditional distribution gives the PBR
                RNS(a,b,x,y) = FreqEst(a,b,x,y)/Q_opt(a,b,x,y)*Pxy(x,y);
                betaNS(a,b,x,y) = RNS(a,b,x,y)*Pxy(x,y);
                
            end
        end
    end
end

% Ideally, the KL divergence minimization should return the global unique 
% minimizer. However, numerical imprecision in the optimization may result
% in an approximate global minimizer, thus resulting in Bell value that is
% not exactly 1 for the Bell-like inequality defined by betaNS. Here, we
% compute using a linear program the maximal Bell value for the
% corresponding (one-way) NS set

[CorrectedBoundNS, Q_opt, DiagMaxBell] = MaxBellValue_NSCond(betaNS, Dir);

% Correction of the PBR due to numerical imprecision in the KL divergence
% minimization

RNS = RNS/CorrectedBoundNS;

logTNS = 0;

% Computation of the log of the test statistics
for x=1:nX
    for y=1:nY
        for a=1:nA
            for b=1:nB
                
                if NABXYTest(a,b,x,y)~=0
                    
                    logTNS = logTNS+NABXYTest(a,b,x,y)*log(RNS(a,b,x,y));
                    
                end
                
            end
        end
    end
end


% Converting the test statistics to the p-value upper bound
pValNS = exp(-logTNS);

% Specification of output filename
ofilename = strcat('AnalyzeBellTestData.mat');

fprintf('\n The p-value upper bound is %.3e\n', pValNS);

if pValNS < pCutoff
    fprintf(' For the specified p-value cutoff of %.3e,', pCutoff);
    if strcmp(Dir,'NS')
        fprintf('\n the hypothesis of two-way NS should be rejected\n\n')
    else
        fprintf('\n the hypothesis of %s one-way NS should be rejected\n\n', Dir);
    end
end

save(ofilename);

