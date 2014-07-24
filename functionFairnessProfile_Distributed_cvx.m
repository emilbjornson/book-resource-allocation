function [finalRates,WBestBeamforming,allPoints] = functionFairnessProfile_Distributed_cvx(H,D,Qsqrt,q,epsilon,lowerPoint,direction,xi,tau,maxIterations)
%Solves the fairness-profile optimization (FPO) problem in Example 2.8
%using the distributed Algorithm 5. The FPO problem searches on a line
%from LOWERPOINT in the direction specified by DIRECTION. The algorithm
%tries to find the intersection between the line and the Pareto boundary of
%the rate region.
%
%The FPO problem is equivalent to finding (g_1,...,) that solves
%
%maximize min_k ( g_k - lowerPoint(k) ) / direction(k) 
%
%subject to     g_k >= lowerPoint(k) for all users k,
%               Power constraints.
%
%The distributed algorithm moves step-by-step along the line using the
%fixed step size TAU. The feasibility of a point on the line is evaluated
%by solving the feasibility problem (2.29) in a distributed manner. Dual
%decomposition is applied to divide the original problem (2.29) into the Kr
%distributed subproblems in (4.30) and a centralized master problem (4.31).
%The subproblems are convex optimization problems while the master problem
%is solved by moving along the subgradient with a step size of XI. 
%
%This is version 1.0.
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%INPUT:
%H             = Kr x Kt*Nt matrix with row index for receiver and column
%                index transmit antennas
%D             = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if
%                j:th transmit antenna can transmit to user k and zero otherwise
%Qsqrt         = N x N x L matrix with matrix-square roots of the L weighting
%                matrices for the L power constraints
%q             = Limits of the L power constraints
%epsilon       = Accuracy of consistency constraints in the subproblems
%lowerPoint    = Start point of the line (must be inside the rate region)
%direction     = Kr x 1 vector with line direction from the start point
%xi            = 
%tau           = Step size along the line (improvement in sum rate at the
%                new point)
%maxIterations = Number of times that the distributed subproblems are
%                solved.
%
%OUTPUT:
%finalRates       = Kr x 1 vector with the best feasible point found by the
%                   algorithm
%WBestBeamforming = Kt*Nt x Kr matrix with beamforming that achieves the
%                   feasible point finalRates
%allPoints        = Kr x maxIterations matrix with the feasible points that
%                   are found by the algorithm. Columns are only non-zero
%                   for iterations when a new better feasible point was
%                   found.


Kr = size(H,1); %Number of users
N = size(H,2); %Number of transmit antennas (in total)
L = size(Qsqrt,3); %Number of power constraints


%Pre-allocation of matrix for storing optimal beamforming
WBestBeamforming = zeros(N,Kr);

%Pre-alloaction of matrix to store the point achieved at each iteration
allPoints = zeros(Kr,maxIterations);

%Pre-allocation of matrices for Lagrange multipliers
Y = zeros(Kr,Kr);
z = zeros(L,1);


%Initialize the step size of the subgradient method for the master problem
xi_current = xi;


%Initialize the QoS requirement that should be solved
candidatePoint = lowerPoint+tau*direction;


%Iterate over the number subproblems that should be solved.
for m = 1:maxIterations
    
    gammavar = 2.^(candidatePoint)-1; %Current QoS requirements
        
    %Preallocation of matrices for storing current values on coupling
    %variables and beamforming in the subproblems.
    Theta_current = zeros(Kr,Kr);
    ThetaTilde_current = zeros(Kr,Kr);
    qk_current = zeros(L,Kr);
    W = zeros(N,Kr);
    
    %Iteration over users (for solving their subproblems)
    for k = 1:Kr
        
        %Solve the subproblem (4.30) for User k using CVX
        cvx_begin
        cvx_quiet(true); % this suppresses screen output from the solver
        
        variable v(N) complex;  %Variable for N x 1 beamforming vector
        variable Thetak(Kr) %Coupling variables Theta_ki for i=1...Kr
        variable ThetaTildek(Kr) %Coupling variables tilde{Theta}_ik for i=1...Kr
        variable qk(L) %Coupling variables q_lk for l=1...L
        
        minimize Y(k,:)*Thetak - transpose(Y(:,k))*ThetaTildek + transpose(z)*qk %Minimize the partial Lagrangian for fixed A and b
        
        %SOCP formulation for the SINR constraint in the k:th subproblem
        hkD = H(k,:)*D(:,:,k);
        imag(hkD*v) == 0;
        real(hkD*v) >= sqrt(gammavar(k))*norm([1; ThetaTildek([1:k-1 k+1:Kr])]);
        
        %Power constraints in the k:th subproblem
        for l = 1:L
            norm(Qsqrt(:,:,l)*v) <= sqrt(qk(l));
            qk(l) <= q(l);
        end
        
        %Interference constraints in the k:th subproblem
        for i=[1:k-1 k+1:Kr]
            norm(H(i, :)*D(:,:,k)*v) <= Thetak(i);
        end
        
        cvx_end
        
        %Store the solution of each subproblem.
        if strfind(cvx_status,'Solved')
            Theta_current(k,:) = Thetak';
            ThetaTilde_current(:,k) = ThetaTildek;
            qk_current(:,k) = qk;
            W(:,k) = v;
        end
        
    end
    
    
    %The diagonal in the coupling variables should be zero, since users are
    %not causing interference to themselves.
    Theta_current = Theta_current-diag(diag(Theta_current));
    ThetaTilde_current = ThetaTilde_current-diag(diag(ThetaTilde_current));

    c = max(sum(qk_current,2)./q); %Calculate scaling factor as in Step 6.
    
    %Scale variables as described in Step 7
    W = W/sqrt(c);
    Theta_current = Theta_current/sqrt(c);
    ThetaTilde_current = ThetaTilde_current/sqrt(c);
    qk_current = qk_current/c;
    

    %Move towards the solution of the master dual problem in (4.31) for the
    %current values on the Lagrange multipliers. This is achieved by the
    %subgradient method described in (4.32) and (4.33).
    Y = Y - xi_current * (ThetaTilde_current-Theta_current); %Update Lagrange multipliers for coupling variables related to SINR constraints
    Y = Y - diag(diag(Y)); %The diagonal is not used
    z = z - xi_current * (q - sum(qk_current,2)); %Update Lagrange multipliers for coupling variables related to power constraints
    Y(Y<0) = 0; %Lagrange multipliers are positive (or non-negative).
    z(z<0) = 0; %Lagrange multipliers are positive (or non-negative).
    
    
    %Check consistency constraints as in Step 9
    interference_consistency = Theta_current-ThetaTilde_current <= epsilon;
    power_consistency = c<1;
    
    %Compute the achieved SINRs as in Step 12.
    SINRs = zeros(Kr,1);
    for k = 1:Kr
        SINRs(k) = norm(H(k,:)*D(:,:,k)*W(:,k))^2/(1+sum(Theta_current(:,k).^2)); 
    end

    
    %Check if the current QoS requirements are satisfied, either since the
    %consistency constraints are satisfied (to accuracy epsilon) or since
    %the achieved SINRs satisfy them.
    if (min(interference_consistency(:)>0) && power_consistency) || min(SINRs>=gammavar)>0
        WBestBeamforming = W; %Store current beamforming
        allPoints(:,m) = log2(1+SINRs); %Store the current operating point
        
        lowerPoint = min(log2(1+SINRs)./direction) * direction; %Update operating point by projecting the current point onto the search curve
        candidatePoint = lowerPoint + tau*direction; %Update the candidate point
        
        xi_current = xi; %Set the step size in the subgradient method to its initial value
        
    else
        %If the QoS requirements are not satisfied: Reduce the step size of
        %the subgradient method and continue with a new iteration.
        factor = (xi/xi_current)^2+1;
        xi_current = xi/sqrt(factor);
        
    end    

end

%Calculate user rates with best feasible beamforming that has been found
channelGains = abs(H*WBestBeamforming).^2;
signalGains = diag(channelGains);
interferenceGains = sum(channelGains,2)-signalGains;
finalRates = log2(1+signalGains./(1+interferenceGains));
