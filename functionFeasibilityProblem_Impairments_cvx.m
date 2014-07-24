function [feasible,Wsolution] = functionFeasibilityProblem_Impairments_cvx(H,D,Qsqrt,q,gammavar,EVMs)
%Solves the feasibility problem with quality-of-service (QoS) constraints
%and transceiver impairments in (4.46). The implementation is based on the
%convex formulation in Corollary 4.3, but has been simplified to only
%handle linear distortion functions.

%For numerical reasons, the implementation considers power minimization
%under the QoS requirements:
%
%minimize   betavar
%subject to SINR_k >= gammavar(k) for all users k,
%           Impairment constraints for all antennas,
%           Power constraints scaled by betavar.
%
%If this optimization problem is feasible and betavar<=1, then the 
%feasibility problem with QoS constraints is also feasible.
%
%This optimization problem is convex. The computational complexity is
%therefore polynomial in the number of users, antennas, and power
%constraints. The implementation can, at least, handle 30 users, 50
%antennas, and 50 power constraints.
%
%This is version 1.0.
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%INPUT:
%H          = Kr x Kt*Nt matrix with row index for receiver and column
%             index transmit antennas
%D          = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if
%             j:th transmit antenna can transmit to user k and zero otherwise
%Qsqrt      = N x N x L matrix with matrix-square roots of the L weighting 
%             matrices for the L power constraints
%q          = Limits of the L power constraints
%gammavar   = Kr x 1 vector with SINR constraints for all users.
%EVMs       = 2 x 1 vector with EVM at the transmit and receive antennas,
%             respectively. These are used along with linear distortion
%             functions.
%
%OUTPUT:
%feasible  = This variable is feasible=true if the feasibility problem is
%            feasible. Otherwise we have feasible=false.
%Wsolution = Kt*Nt x Kr matrix with beamforming achieved by the power
%            minimization problem. This matrix is empty if this problem is
%            infeasible.


Kr = size(H,1); %Number of users
N = size(H,2); %Number of transmit antennas (in total)
L = size(Qsqrt,3); %Number of power constraints

EVMt = EVMs(1);
EVMr = EVMs(2);
EVMrPlusNoise = sqrt(1+EVMr^2);


%Solve the power minimization under QoS requirements problem using CVX
cvx_begin
cvx_quiet(true); % This suppresses screen output from the solver

variable W(N,Kr) complex;  %Variable for N x Kr beamforming matrix
variable t(N); %Auxiliary variable for N x 1 distortion magnitude at the transmit antennas
variable betavar %Scaling parameter for power constraints

minimize betavar %Minimize the power indirectly by scaling power constraints

subject to

%SINR constraints (Kr constraints)
for k = 1:Kr
    
    %Channels of the signal intended for user i when it reaches user k
    hkD = zeros(Kr,N);
    for i = 1:Kr
        hkD(i,:) = H(k, :)*D(:,:,i);
    end
    
    imag(hkD(k, :)*W(:, k)) == 0; %Useful link is assumed to be real-valued
    
    %SOCP formulation for SINR constraint of user k - Adapted for
    %transceiver impairments at the transmitter and receiver. Note that the
    %impact of impairments at the receiver can be included directly in the
    %SOCP when the linear distortion functions are considered.
    real(hkD(k, :)*W(:, k)) >= sqrt(gammavar(k))*norm([1; EVMrPlusNoise*diag(hkD([1:k-1 k+1:Kr],:)*W(:,[1:k-1 k+1:Kr])); EVMr*hkD(k,:)*W(:,k); abs(hkD(k,:)') .* t ]);

end

%The distortion magnitude at each transmit antenna should be smaller or
%equal to the auxiliary distortion magnitude variable t. Equality will be
%achieved at the optimal solution.
for n=1:N
    EVMt*norm(W(n,:)) <= t(n);
    t(n) >= 0;
end

%Power constraints (L constraints) scaled by the variable betavar
for l = 1:L
    norm([Qsqrt(:,:,l)*W diag(Qsqrt(:,:,l)).*t],'fro') <= betavar*sqrt(q(l));
end

betavar >= 0; %Power constraints must be positive

cvx_end



%Analyze result and prepare the output variables.
if isempty(strfind(cvx_status,'Solved')) %Both power minimization problem and feasibility problem are infeasible.
    feasible = false;
    Wsolution = [];
elseif betavar>1 %Only power minimization problem is feasible.
    feasible = false;
    Wsolution = W;
else %Both power minimization problem and feasibility problem are feasible.
    feasible = true;
    Wsolution = W;
end
