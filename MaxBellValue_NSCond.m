% MaxBellValue_NSCond
% requires: yalmip, mosek

% author: Yeong-Cherng Liang
%
%    [BellVal, Q_opt, diagnostics] = MaxBellValue_NSCond(beta, Dir)
%    maximizes the Bell value of a bipartite Bell-like inequality specified 
%    by beta(a,b,x,y) with respect to correlations from the one-way or 
%    two-way no-signaling (NS) set. For the one-way NS set, the direction 
%    of NS constraints are to be specified by "Dir"; if only two arguments 
%    are specified, the usual NS set is assumed. The maximal Bell-value for 
%    the speicifed set is returned as "BellVal" alongside the maximizing
%    conditional distribution Q_opt, and the diagnostics information from
%    the optimization "diagonistics"

%% Copyright (C) 2025 Yeong-Cherng Liang, last modified on 17 Mar 2025

function [BellVal, Q_opt, diagnostics] = MaxBellValue_NSCond(beta, Dir)

    % Validate input
    if ndims(beta) ~= 4
        error('beta(a, b, x, y) must be a 4D array.');
    end
    [nA, nB, nX, nY] = size(beta);

    % Flatten beta to make it a column vector for the convenience of 
    % subsequent computation
    beta_flat = beta(:);

    % Define the sdpvar corresponding to the one-way or two-say 
    % no-signaling correlation, i.e., the conditional distribution Q_flat
    Q_flat = sdpvar(size(beta_flat,1), 1);


    % Objective: Maximize the Bell value 
    objective = sum(beta_flat.*Q_flat);

    Q = reshape(Q_flat, [2 2 2 2]);

    % Constraints: Q is a valid probability distribution
    constraints = [
        % All entries of Q must be non-negative
        Q_flat >= 0,              
        % Imposing the normalization constraint for Q_flat
        squeeze(sum(Q,[1 2]))==1, 
    ];

    % Imposing one-way no-singaling constraints from B to A
    if nargin<3 || (nargin==3 && strcmp(Dir,"B2A"))
        for a = 1:nA
            for x = 1:nX
                for y = 1:nY
                    for y_prime = 1:nY
                        if y_prime ~= y
                            constraints = [constraints, sum(Q(a,:,x,y)-Q(a,:,x,y_prime))==0];
                        end
                    end
                end
            end
        end
    end

    % Imposing one-way no-singaling constraints from A to B
    if nargin<3 || (nargin==3 && strcmp(Dir,"A2B"))
        for b = 1:nB
            for y = 1:nY
                for x = 1:nX
                    for x_prime = 1:nX
                        if x_prime ~= x
                            constraints = [constraints, sum(Q(:,b,x,y)-Q(:,b,x_prime,y))==0];
                        end
                    end
                end
            end
        end
    end

    % Setting optimization parameters (verbose set to 0 to run the
    % optimization quitely)
    options = sdpsettings('solver', 'mosek', ...
    'mosek.MSK_DPAR_INTPNT_TOL_PFEAS',1e-9, ...
    'mosek.MSK_DPAR_INTPNT_TOL_DFEAS',1e-9, ...
    'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP',1e-9, ...
    'verbose', 0);
    
    % We maximize the objective function to minimize KL divergence
    diagnostics = optimize(constraints, -objective, options);

    % Check for solver success
    if diagnostics.problem == 0

        % Extract the solution
        Q_opt_flat = value(Q_flat);
        Q_opt = reshape(Q_opt_flat, [nA, nB, nX, nY]);
        BellVal = value(objective);

    else
        fprintf('Optimization failed: %s\n', diagnostics.info);
    end
    yalmip('clear');
end