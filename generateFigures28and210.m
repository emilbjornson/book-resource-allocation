%This m-file generates Figure 2.8 and Figure 2.10 in "Optimal Resource
%Allocation in Coordinated Multi-Cell Systems" by Emil Björnson and Eduard
%Jorswieck.
%
%This is version 1.2. (Last edited: 2014-07-23)
%
%Figure 2.8 illustrates the Polyblock Outer Approximation (PA) algorithm. 
%The sum information rate is maximized by approximating the performance
%region (around the optimal point) from above using a polyblock. The
%approximation is iteratively refined and parts that cannot contain the
%optimal point are removed.
%
%Figure 2.10 illustrates the Branch-Reduce-and-Bound (BRB) algorithm.
%The sum information rate is maximized by approximating the Pareto boundary
%of the performance region (around the optimal point) using a set of
%disjoint boxes. The approximation is iteratively refined and parts that
%cannot contain the optimal point are removed.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the book.

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
pathlossOverNoise = 1/2*eye(Kr)+1/2*ones(Kr,Kr);
H = sqrt(kron(pathlossOverNoise,ones(1,Nt))).*(randn(Kr,Nt*Kt)+1i*randn(Kr,Nt*Kt))/sqrt(2);



%%Part 1: Calculate sample points on the Pareto boundary using Approach 1
%%in Section 3.3. The samples are computed by applying Theorem 3.4 on fine
%%grid of weighting vectors.

deltaFPO = 0.1; %Accuracy of the line searches in the FPO algorithm

nbrOfSamplesFPO = 101; %Number of sample points
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
        disp(['Progress: ' num2str(m) ' out of ' num2str(nbrOfSamplesFPO) ' search directions']); 
    end
    
    %Search on a line from the origin in each of the equally spaced search
    %directions. The upper point is either on Pareto boundary (unlikely) or
    %outside of the rate region.
    lowerPoint = zeros(Kr,1);
    upperPoint = sum(utopia) * profileDirections(:,m);
    
    %Find the intersection between the line and the Pareto boundary by
    %solving an FPO problem.
    finalInterval = functionFairnessProfile_cvx(H,D,Qsqrt,q,deltaFPO,lowerPoint,upperPoint);
    
    ParetoBoundary(:,m) = finalInterval(:,1); %Save the best feasible point found by the algorithm
end

%Search for the sample points that maximize the sum rate. Note that the
%accuracy is improved when the number of sample points is increased.
[~,sumrateIndex] = max(ParetoBoundary(1,:)+ParetoBoundary(2,:));



%%Part 2: Run the PA and BRB algorithms to search for the point that
%%maximizes the sum rate.

disp('Progress: Running PA and BRB algorithms'); %Output of the simulation progress

%Define optimization problem as maximizing the unweighted sum rate.
problemMode = 1;
weights = ones(Kr,1);

saveBoxes = 1; %Tell PA and BRB algorithms to return the set of vertices/boxes, to enable plotting.

epsilon = 0.001; %Accuracy of the optimal sum rate in the BRB algorithm
deltaBRB = 0.001; %Accuracy of the line searches in the BRB algorithm
deltaPA=0.001; %Accuracy of the line searches in the PA algorithm

%Since the purpose of this figure is to illustrate the algorithms, we will
%not run the algorithms to convergence but only 40 iterations.
maxIterations = 30; %Maximal number of iterations of the algorithm
maxFuncEvaluations = 500; %Maximal number of convex problems solved in the algorithm

%Define the origin, which is used as the best feasible point known in
%advance (in both algorithms) and as the lower corner of the initial box in
%the BRB algorithm.
origin = zeros(Kr,1);

%Run the PA algorithm
[bestFeasible,~,~,~,vertices] = functionPAalgorithm_cvx(H,D,Qsqrt,q,utopia,weights,deltaPA,epsilon,maxIterations,maxFuncEvaluations,origin,problemMode,saveBoxes);

%Run the BRB algorithm
[bestFeasibleBRB,~,~,~,boxes] = functionBRBalgorithm_cvx(H,D,Qsqrt,q,origin,utopia,weights,deltaBRB,epsilon,maxIterations,maxFuncEvaluations,origin,problemMode,saveBoxes);

disp('Progress: Finalizing'); %Output of the simulation progress



%%Plot the progress of PA and BRB algorithms.

whichIterations=[0 1 2 10 20 30]; %Which 6 iterations will be shown.

%PA algorithm
figure(1);

for m = 1:6
    %Plot the rate region and point for maximal sum rate.
    subplot(3,2,m); hold on; box on;
    plot(ParetoBoundary(1,:),ParetoBoundary(2,:),'k');
    plot(ParetoBoundary(1,sumrateIndex),ParetoBoundary(2,sumrateIndex),'k*');
    
    if whichIterations(m) == 0
        %First subplot shows the initial box
        plot([origin(1) origin(1) utopia(1) utopia(1) origin(1)],[origin(2) utopia(2) utopia(2) origin(2) origin(2)],'b','LineWidth',0.75);
        
        legend('Pareto Boundary','Maximal Sum Utility','Location','SouthWest');
        xlabel('Initialization (PA algorithm)');
    else
        %m:th subplot shows the polyblock at iteration whichIterations(m).
        %We first need to convert the set of vertices into a polyblock.
        vertexSet = vertices{whichIterations(m)}.set;
        vertexSetSorted = [origin(1) max(vertexSet(2,:)); sortrows(vertexSet',1); max(vertexSet(1,:)) origin(2)];
        polyblockCurve = zeros(2*size(vertexSetSorted,1),Kr);
        for n = 1:size(vertexSetSorted,1)-1
            polyblockCurve(2*(n-1)+1,:) = vertexSetSorted(n,:);
            polyblockCurve(2*(n-1)+2,:) = [vertexSetSorted(n,1) vertexSetSorted(n+1,2)];
        end
        polyblockCurve(end-1,:) = vertexSetSorted(end,:);
        
        plot(polyblockCurve(:,1),polyblockCurve(:,2),'b','LineWidth',0.75); %Plot the polyblock
        
        xlabel(['Iteration ' num2str(whichIterations(m))]);
    end
end

%BRB algorithm
figure(2);

for m = 1:6
        %Plot the rate region and point for maximal sum rate.
    subplot(3,2,m); hold on; box on;
    plot(ParetoBoundary(1,:),ParetoBoundary(2,:),'k');
    plot(ParetoBoundary(1,sumrateIndex),ParetoBoundary(2,sumrateIndex),'k*');
    
    if whichIterations(m) == 0
        %First subplot shows the initial box
        plot([origin(1) origin(1) utopia(1) utopia(1) origin(1)],[origin(2) utopia(2) utopia(2) origin(2) origin(2)],'b','LineWidth',0.75);
        
        legend('Pareto Boundary','Maximal Sum Utility','Location','SouthWest');
        xlabel('Initialization (BRB algorithm)');
    else
        %m:th subplot shows the set of boxes at iteration whichIterations(m).
        boxesLowerCorners = boxes{whichIterations(m)}.lowerCorners;
        boxesUpperCorners = boxes{whichIterations(m)}.upperCorners;
        
        for boxInd = 1:size(boxesUpperCorners,2)
            plot([boxesLowerCorners(1,boxInd) boxesLowerCorners(1,boxInd) boxesUpperCorners(1,boxInd) boxesUpperCorners(1,boxInd) boxesLowerCorners(1,boxInd)],[boxesLowerCorners(2,boxInd) boxesUpperCorners(2,boxInd) boxesUpperCorners(2,boxInd) boxesLowerCorners(2,boxInd) boxesLowerCorners(2,boxInd)],'b','LineWidth',0.75);
        end
        
        xlabel(['Iteration ' num2str(whichIterations(m))]);
    end
end
