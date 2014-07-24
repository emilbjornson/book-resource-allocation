%This m-file generates Figure 4.5 in "Optimal Resource Allocation in
%Coordinated Multi-Cell Systems" by Emil Björnson and Eduard Jorswieck.
%
%This is version 1.2. (Last edited: 2014-07-23)
%
%Figure 4.5 illustrates the convergence of the distributed resource
%allocation algorithm in Algorithm 5. Max-min fairness is the system
%utility function and large percentage of the optimal performance is
%achieved after a few iterations.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the book. The code can be improved by
%having adaptive step sizes in the algorithm for distributed
%fairness-profile optimization; this will improve the convergence for
%certain channel realizations.

close all;
clear all;


%%Simulation parameters

Kt = 2; %Number of base stations (BSs)
Kr = 2; %Number of users (in total)
Nt = 2; %Number of antennas per BS

rng('shuffle'); %Initiate the random number generators with a random seed

%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));
%rand('state',sum(100*clock));

PdB = 10; %SNR (in dB)
P = 10.^(PdB/10); %SNR

%D-matrices for global joint transmission
D = repmat(eye(Kt*Nt),[1 1 Kr]);

%Definition of per base station power constraints
L = Kt;
Q = zeros(Kt*Nt,Kt*Nt,L); %The L weighting matrices
Qsqrt = zeros(Kt*Nt,Kt*Nt,L); %The matrix-square root of the L weighting matrices
for l = 1:L
    Q((l-1)*Nt+1:l*Nt,(l-1)*Nt+1:l*Nt,l) = (1/P)*eye(Nt);
    Qsqrt((l-1)*Nt+1:l*Nt,(l-1)*Nt+1:l*Nt,l) = sqrt(1/P)*eye(Nt);
end
q = ones(L,1); %Limits of the L power constraints. Note that these have been normalized.

%Generation of normalized Rayleigh fading channel realizations with unit
%variance to closest base station and 1/3 as variance from other base station.
pathlossOverNoise = 2/3*eye(Kr)+1/3*ones(Kr,Kr);
H = sqrt(kron(pathlossOverNoise,ones(1,Nt))).*(randn(Kr,Nt*Kt)+1i*randn(Kr,Nt*Kt))/sqrt(2);



%%Part 1: Run the distributed Algorithm 5 to solve the fairness-profile
%%optimization (FPO) problem.

disp('Progress: Running distributed algorithm'); %Output of the simulation progress

%Problem formulation
lowerPoint = zeros(Kr,1); %Start the search in the origin
direction = ones(Kr,1)/Kr; %Search in the direction of equal user performance (unweighted max-min fairness)

%Number of iterations in the algorithm in terms of number of times that the
%distributed subproblems are solved. This number is small for illustrative
%purposes.
maxIterations = 10;

%Parameters in the distributed algorithm
tau = 0.1; %Step size along the search curve
xi = 5; %Step size of the subgradient method in the master problem
epsilon = 0.01; %Accuracy of consistency constraints in the subproblems

%Run the distributed algorithm
[finalRates,WBestBeamforming,allPoints] = functionFairnessProfile_Distributed_cvx(H,D,Qsqrt,q,epsilon,lowerPoint,direction,xi,tau,maxIterations);



%%Part 2: The achievable points from the distributed algorithm will be
%%illustrated and compared towards to Pareto boundary of the rate region.
%%Sample points on the Pareto boundary are computed by applying
%%Theorem 3.4 on fine grid of weighting vectors.

delta = 0.1; %Accuracy of the line searches in the FPO algorithm

nbrOfSamplesFPO = 1001; %Number of sample points
ParetoBoundary = zeros(Kr,nbrOfSamplesFPO); %Pre-allocation of sample matrix


%Computation of the utopia point using MRT with full transmit power,
%which is optimal when there is only one active user.

wMRT = functionMRT(H,D); %Compute normalized beamforming vectors for MRT

utopia = zeros(Kr,1);
for k = 1:Kr
    W = sqrt(P)*[wMRT(1:2,k)/norm(wMRT(1:2,k)); wMRT(3:4,k)/norm(wMRT(3:4,k))]; %Scale to use all available power
    utopia(k) = log2(1+abs(H(k,:)*W)^2);
end


%Generate a grid of equally spaced search directions
interval = linspace(0,1,nbrOfSamplesFPO);
profileDirections = [interval; 1-interval];

for m = 1:nbrOfSamplesFPO

    %Output of the simulation progress - since this part takes a long time
    if mod(m,10) == 1
        disp(['Progress: Generating boundary, ' num2str(m) ' out of ' num2str(nbrOfSamplesFPO) ' search directions']); 
    end
    
    %Search on a line from the origin in each of the equally spaced search
    %directions. The upper point is either on Pareto boundary (unlikely) or
    %outside of the rate region.
    lowerPoint = zeros(Kr,1);
    upperPoint = sum(utopia) * profileDirections(:,m);
    
    %Find the intersection between the line and the Pareto boundary by
    %solving an FPO problem.
    finalInterval = functionFairnessProfile_cvx(H,D,Qsqrt,q,delta,lowerPoint,upperPoint);
    
    ParetoBoundary(:,m) = finalInterval(:,1); %Save the best feasible point found by the algorithm
end

%Search for the sample points that maximize the max-min fairness and sum
%rate, respectively.
[~,maxminIndex] = max(min(ParetoBoundary,[],1));
[~,sumrateIndex] = max(ParetoBoundary(1,:)+ParetoBoundary(2,:));




%Plot the rate region and the result of the distributed algorithm
figure; hold on; box on;

%Plot sample points that maximize sum rate and max-min fairness, respectively.
plot(ParetoBoundary(1,sumrateIndex),ParetoBoundary(2,sumrateIndex),'ro');
plot(ParetoBoundary(1,maxminIndex),ParetoBoundary(2,maxminIndex),'k*');

plot(ParetoBoundary(1,:),ParetoBoundary(2,:),'k-','LineWidth',1); %Plot sample points of Pareto boundary

%Plot search curve
upperPoint = min(utopia./direction) * direction;
plot([0 upperPoint(1)],[0 upperPoint(2)],'k--');

%Plot points achieved by the distributed algorithm
for m = 1:maxIterations
    
    %Check at which iterations a new better feasible point is found
    if sum(allPoints(1,m))>0
        plot(allPoints(1,m),allPoints(2,m),'k+');  %Plot the new feasible point
        
        %Calculate appropriate position of text
        xTextLocation = max([allPoints(1,m)-0.25 0.1]);
        
        if allPoints(2,m) > 0.4
            yTextLocation=allPoints(2,m)-0.25;
        else
            yTextLocation=allPoints(2,m)+0.25;
        end
        
        text(xTextLocation,yTextLocation,['Iteration ' num2str(m)]); %Plot text with iteration number
    end
end

xlabel('log_2(1+SINR_1) [bits/channel use]')
ylabel('log_2(1+SINR_2) [bits/channel use]')

legend('Max Arithmetic Mean','Max-Min Fairness','Location','NorthEast')
