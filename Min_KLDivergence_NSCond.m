% Min_KLDivergence_NSCond
% requires: yalmip, mosek

% author: Yeong-Cherng Liang
%
%    [Q_opt, diagnostics] = Min_KLDivergence_NSCond(P, Pxy, Dir)
%    minimizes the KL divergence of a given P(a,b|x,y) with priors Pxy
%    to the NS set. The minimizer is returned as the unconditional 
%    Q_opt(a,b,x,y). The NS set is enforced here by imposing all the NS 
%    constraints in the optimization problem.
%

%% Copyright (C) 2025 Yeong-Cherng Liang, last modified on 17 Mar 2025

function [Q_opt, diagnostics] = Min_KLDivergence_NSCond(P, Pxy, Dir)

    % Validate input
    if ndims(P) ~= 4
        error('P(a, b, x, y) must be a 4D array.');
    end
    [nA, nB, nX, nY] = size(P);

    if isscalar(Pxy)
        % Compute the joint distribution for a flat input distribution 
        Pjoint = Pxy*P;
    else
        % Compute the joint distribution for a general input distribution 
        Pjoint = zeros(nA,nB,nX,nY);
        
        for a = 1:nA
            for b = 1:nB
                for x = 1:nX
                    for y = 1:nY
                        
                        Pjoint(a,b,x,y) = Pxy(x,y)*P(a,b,x,y);
                        
                    end
                end
            end
        end
    end

    
    % Flatten P to make it a column vector for the convenience of 
    % subsequent computation
    P_flat = Pjoint(:);

    % Define the sdpvar corresponding to the minimizing joint 
    % (unconditional) distribution Q_flat and the auxiliary variables u_flat
    Q_flat = sdpvar(size(P_flat,1), 1);
    u_flat = sdpvar(size(P_flat,1), 1);


    % Objective: Minimize KL divergence 
    %     KL_div = sum(P_flat.*(log(P_flat)-log(Q_flat))),
    % which we write in its conic form after dropping the constant term
    % and introducing auxiliary variables u

    objective = sum(P_flat.*u_flat);
    
    % Constraints: Q is a valid probability distribution
    constraints = [
        % All entries of Q must be non-negative
        Q_flat >= 0,              
        % Imposing the normalization constraint for Q_flat
        sum(Q_flat)==1, 
        % Imposing the relations between Q_flat and the auxiliary variables
        Q_flat >= exp(u_flat),
    ];

    Q = reshape(Q_flat, [2 2 2 2]);

    % Imposing one-way no-singaling constraints from B to A
    if nargin<3 || (nargin==3 && strcmp(Dir,"B2A"))
        for a = 1:nA
            for x = 1:nX
                for y = 1:nY
                    for y_prime = 1:nY
                        if y_prime ~= y
                            if isscalar(Pxy)
                                constraints = [constraints, sum(Q(a,:,x,y)-Q(a,:,x,y_prime))==0];
                            else
                                constraints = [constraints, sum(Q(a,:,x,y)/Pxy(x,y)-Q(a,:,x,y_prime)/Pxy(x,y_prime))==0];
                            end
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
                            if isscalar(Pxy)
                                constraints = [constraints, sum(Q(:,b,x,y)-Q(:,b,x_prime,y))==0];
                            else
                                constraints = [constraints, sum(Q(:,b,x,y)/Pxy(x,y)-Q(:,b,x_prime,y)/Pxy(x_prime,y))==0];
                            end
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

    else
        fprintf('Optimization failed: %s\n', diagnostics.info);
    end
    yalmip('clear');
end