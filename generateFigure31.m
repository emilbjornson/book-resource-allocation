%This m-file generates Figure 3.1 in "Optimal Resource Allocation in
%Coordinated Multi-Cell Systems" by Emil Björnson and Eduard Jorswieck.
%
%This is version 1.2. (Last edited: 2014-07-23)
%
%Figure 3.1 illustrates the performance regions for random channel
%realizations in a scenario with two users and global joint transmission.
%The Pareto boundaries are approximated in two ways: (1) Solving FPO over a
%grid of weighting vectors; and (2) Generating sample points using a
%beamforming parametrization.
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
pathlossOverNoise = 2/3*eye(Kr)+1/3*ones(Kr,Kr);
H = sqrt(kron(pathlossOverNoise,ones(1,Nt))).*(randn(Kr,Nt*Kt)+1i*randn(Kr,Nt*Kt))/sqrt(2);




%%Approach 1: Calculate sample points on the Pareto boundary by applying
%%Theorem 3.4 on fine grid of weighting vectors.

disp('Progress: Running Approach 1'); %Output of the simulation progress

delta = 0.1; %Accuracy of the line searches in the FPO algorithm

nbrOfSamplesFPO = 1001; %Number of sample points
ParetoBoundary = zeros(Kr,nbrOfSamplesFPO); %Pre-allocation of sample matrix


%Computation of the utopia point using MRT with full transmit power, 
%which is optimal when there is only one active user. Note that the utopia
%point is only calculated in this way under global joint transmission.

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
        disp(['Progress: Running Approach 1, ' num2str(m) ' out of ' num2str(nbrOfSamplesFPO) ' search directions']); 
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

%Add the single-user points to the Pareto boundary.
ParetoBoundary = [[0; utopia(2)] ParetoBoundary [utopia(1); 0]];


%Search for the sample points that maximize the sum rate, geometric mean,
%harmonic mean, and max-min fairness, respectively. Note that the accuracy
%of these points improves when the number of sample points is increased.
[~,sumrateIndex] = max(ParetoBoundary(1,:)+ParetoBoundary(2,:));
[~,geomeanIndex] = max(ParetoBoundary(1,:).*ParetoBoundary(2,:));
[~,harmmeanIndex] = max(harmmean(ParetoBoundary,1));
[~,maxminIndex] = max(min(ParetoBoundary,[],1));



%%Approach 2: Calculate sample points of the whole performance region using
%%the beamforming parametrization in Section 3.2.3.

disp('Progress: Running Approach 2'); %Output of the simulation progress

nbrOfLambdaSamples = 101; %Number of sample points of lambda-parameters
nbrOfMuSamples = 101; %Number of sample points of mu-parameters

lambda1Range = linspace(1e-12,1-1e-12,nbrOfLambdaSamples); %Generate equally spaced samples for lambda1
mu1Range = linspace(1e-12,1-1e-12,nbrOfMuSamples); %Generate equally spaced samples for mu1

%Pre-allocated matrices for storing the sample points for the users
samplePointsUser1 = zeros(length(lambda1Range),length(mu1Range));
samplePointsUser2 = zeros(length(lambda1Range),length(mu1Range));


%Iteration over grid of lambda_k and mu_l values
for lIter = 1:length(lambda1Range)
    
    %Sample values on lambda_1 and lambda_2. These have unit sum, as
    %shown by Corollary 3.2.
    lambda1 = lambda1Range(lIter);
    lambda2 = 1-lambda1;
    
    for mIter = 1:length(mu1Range)
        
        %Sample values on mu_1 and mu_2. These have unit sum, as
        %shown by Corollary 3.2.
        mu1 = mu1Range(mIter);
        mu2 = 1-mu1;
        
        %Compute unnormalized versions of the beamforming directions in (3.16)
        w1 = (mu1*Q(:,:,1)/q(1)+mu2*Q(:,:,2)/q(2)+lambda2*D(:,:,1)*H(2,:)'*H(2,:)*D(:,:,1))\(D(:,:,1)*H(1,:)');
        w2 = (mu1*Q(:,:,1)/q(1)+mu2*Q(:,:,2)/q(2)+lambda1*D(:,:,2)*H(1,:)'*H(1,:)*D(:,:,2))\(D(:,:,2)*H(2,:)');
        
        %Compute gamma_k in (3.19)
        gammas = [real(lambda1*H(1,:)*D(:,:,1)*w1) real(lambda2*H(2,:)*D(:,:,2)*w2)];
        
        %Compute normalized beamforming directions as in (3.16)
        V = [w1/norm(w1) w2/norm(w2)];
        
        %Generate M matrix as in (3.20)
        M = zeros(Kr,Kr);
        for i = 1:Kr
            for j = 1:Kr
                if i == j
                    M(i,i) = abs(H(i,:)*D(:,:,i)*V(:,i)).^2;
                else
                    M(i,j) = -gammas(j)*abs(H(j,:)*D(:,:,i)*V(:,i)).^2;
                end
            end
        end
        
        %Compute power allocation as in (3.17)
        p = gammas/M;
        
        %Compute final beamforming matrix according to Theorem 3.3
        W = [sqrt(p(1))*w1 sqrt(p(2))*w2];
        
        %Scale the beamforming matrix according to Corollary 3.2 to ensure feasibility
        maxPowerAmongBSs = max([norm(Qsqrt(:,:,1)*W,'fro').^2 norm(Qsqrt(:,:,2)*W,'fro').^2]);
        W = W/sqrt(maxPowerAmongBSs);
        
        
        %Calculate and save user rates with the achieved beamforming
        channelGains = abs(H*W).^2;
        samplePointsUser1(lIter,mIter) = log2(1+channelGains(1,1)/(channelGains(1,2)+1));
        samplePointsUser2(lIter,mIter) = log2(1+channelGains(2,2)/(channelGains(2,1)+1));
        
    end
end




%Plot the rate region generated by the two approaches.
figure; hold on; box on;

plot(ParetoBoundary(1,:),ParetoBoundary(2,:),'k-','LineWidth',1); %Plot result of Approach 1
plot(samplePointsUser1(:),samplePointsUser2(:),'.','Color',[.7 .7 .7]); %Plot result of Approach 2

%Plot sample points that maximize sum rate, geometric mean,
%harmonic mean, and max-min fairness, respectively.
plot(ParetoBoundary(1,sumrateIndex),ParetoBoundary(2,sumrateIndex),'ro');
plot(ParetoBoundary(1,geomeanIndex),ParetoBoundary(2,geomeanIndex),'sg');
plot(ParetoBoundary(1,harmmeanIndex),ParetoBoundary(2,harmmeanIndex),'mv');
plot(ParetoBoundary(1,maxminIndex),ParetoBoundary(2,maxminIndex),'k*');

legend('Pareto Boundary based on FPO','Based on Beamforming Parametrization','Max Arithmetic Mean','Max Geometric Mean','Max Harmonic Mean','Max-Min Fairness','Location','SouthWest')

xlabel('log_2(1+SINR_1) [bits/channel use]')
ylabel('log_2(1+SINR_2) [bits/channel use]')
