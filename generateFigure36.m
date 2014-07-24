%This m-file generates Figure 3.6 in "Optimal Resource Allocation in
%Coordinated Multi-Cell Systems" by Emil Björnson and Eduard Jorswieck.
%
%This is version 1.2. (Last edited: 2014-07-23)
%
%Figure 3.6 illustrates the average sum information rate (over channel
%realizations) with SLNR-MAX, MRT, and ZFBF in a 4-user MISO interference
%channel, as a function of the transmit power. MRT and ZFBF are good
%beamforming directions at low and high SNR, respectively, while SLNRMAX
%shows good performance in the entire SNR range.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the book.

close all;
clear all;


%%Simulation parameters

%Should optimal beamforming be calculated? (true /false)
%(Note that it is the calculation of optimal beamforming is the main source
%of computational complexity - the running time varies difference between
%seconds and hours)
computeOptimalBeamforming = true;

Kt = 4; %Number of base stations (BSs)
Kr = 4; %Number of users (in total)
Nt = 4; %Number of antennas per BS

nbrOfRealizations = 10; %Number of channel realizations

rng('shuffle'); %Initiate the random number generators with a random seed

%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));
%rand('state',sum(100*clock));

%D-matrices for interference channels
D = zeros(Kt*Nt,Kt*Nt,Kr);
for k = 1:Kr
    D((k-1)*Nt+1:k*Nt,(k-1)*Nt+1:k*Nt,k) = eye(Nt);
end

%Combined channel matrix will be (Kr x Kt*Nt). This matrix gives the
%normalized variance of each channel element
channelVariances = 1/2+1/2*kron(eye(Kr),ones(1,Nt));

%User weights for (unweighted) sum rate computation
weights = ones(Kr,1);

PdB = -10:1:30; %SNR range (in dB)
P = 10.^(PdB/10); %SNR range



%%System parameters - Special for optimal beamforming (using BRB algorithm)

%Optimal beamforming have its own (sparser) SNR range, since it is the
%Main source of computational complexity
PdB_BRB = -10:5:30; %SNR range (in dB) for optimal beamforming
P_BRB = 10.^(PdB_BRB/10); %SNR range for optimal beamforming

problemMode = 1; %Tells the BRB algorithm to maximize the (weighted) sum rate. (Proportional fairness is given by mode=2)

epsilonBRB = 0.001; %Accuracy of the optimal sum rate in the BRB algorithm
deltaBRB = 0.1; %Accuracy of the line searches in the BRB algorithm
maxIterations = 2000; %Maximal number of iterations of the algorithm
maxFuncEvaluations = 3000; %Maximal number of convex problems solved in the algorithm



%%Pre-allocation of matrices


%Matrices for saving sum rates with different beamforming stategies
sumrateMRT = zeros(length(P),nbrOfRealizations);
sumrateZFBF = zeros(length(P),nbrOfRealizations);
sumrateSLNRMAX = zeros(length(P),nbrOfRealizations);

if computeOptimalBeamforming == true
    sumrateOPTIMAL = zeros(length(P_BRB),nbrOfRealizations);
end

%Pre-generation of Rayleigh fading channel realizations (unit variance)
%(Note that all BSs are assumed to affect all users. If not, the
%corresponding elements of Hall should be set to zero).
Hall = (randn(Kr,Nt*Kt,nbrOfRealizations)+1i*randn(Kr,Nt*Kt,nbrOfRealizations))/sqrt(2);



%Iteration over channel realizations
for m = 1:nbrOfRealizations
    
    %Output of the simulation progress - since this part takes a long time
    disp(['Progress: Initializing channel realization ' num2str(m) ' out of ' num2str(nbrOfRealizations)]); 
    
    %Generate channel matrix for m:th realization
    H = sqrt(channelVariances) .* Hall(:,:,m);
    
    
    %Compute normalized beamforming vectors for MRT
    wMRT = functionMRT(H,D);
    
    %Compute normalized beamforming vectors for ZF
    wZFBF = functionZFBF(H,D);
    
    
    %Pre-allocate matrix for saving good feasible starting points for the
    %BRB algorithm, based on heuristic beamforming
    if computeOptimalBeamforming == true
        bestfeasibleRates = zeros(Kr,length(P_BRB));
        pind_BRB = 1;
    end
    
    
    %Iteration over transmit powers
    for pind = 1:length(P)
        
        %Compute normalized beamforming vectors for SLNR-MAX
        %(Note that it varies with the transmit power)
        wSLNRMAX = functionSLNRMAX(H,P(pind)*ones(Kr,1),D);
        
        
        %Calculate power allocation with MRT (using Theorem 3.5)
        rhos = diag(diag(abs(H*wMRT).^2)); %Interference channel: BS j serves user j
        powerAllocationMRT = functionHeuristicPowerAllocation(rhos,P(pind)*ones(Kt,1),weights); %Interference channel: Allocate full power to each user

        %Calculate sum rate with MRT
        W = kron(sqrt(powerAllocationMRT),ones(Nt,1)).*wMRT;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumrateMRT(pind,m) = weights'*rates;
        
        
        %Calculate power allocation with ZFBF (using Theorem 3.5)
        rhos = diag(diag(abs(H*wZFBF).^2)); %Interference channel: BS j serves user j
        powerAllocationwZFBF = functionHeuristicPowerAllocation(rhos,P(pind)*ones(Kt,1),weights); %Interference channel: Allocate full power to each user
        
        %Calculate sum rate with ZFBF
        W = kron(sqrt(powerAllocationwZFBF),ones(Nt,1)).*wZFBF;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumrateZFBF(pind,m) = weights'*rates;
        
        
        %Calculate power allocation with SLNR-MAX (using Theorem 3.5)
        rhos = diag(diag(abs(H*wSLNRMAX).^2)); %Interference channel: BS j serves user j
        powerAllocationwSLNRMAX = functionHeuristicPowerAllocation(rhos,P(pind)*ones(Kt,1),weights); %Interference channel: Allocate full power to each user
        
        %Calculate sum rate with SLNR-MAX
        W = kron(sqrt(powerAllocationwSLNRMAX),ones(Nt,1)).*wSLNRMAX;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumrateSLNRMAX(pind,m) = weights'*rates;
        
        
        %Save the rates of SLNR-MAX to use as starting point when
        %calculating Optimal beamforming
        if computeOptimalBeamforming == true && P(pind)==P_BRB(pind_BRB)
            bestfeasibleRates(:,pind_BRB) = rates;
            pind_BRB = pind_BRB+1;
        end
        
        
        
    end
    
    if computeOptimalBeamforming == true
        
        for pind = 1:length(P_BRB)
            
    %Output of the simulation progress - since this part takes a long time
    disp(['Progress: BRB Algorithm, ' num2str(pind) ' out of ' num2str(length(P_BRB)) ' transmit powers']); 
    
            
            %Definition of per base station power constraints
            L = Kt;
            Qsqrt = zeros(Kt*Nt,Kt*Nt,L); %The matrix-square root of the L weighting matrices
            for l = 1:L
                Qsqrt((l-1)*Nt+1:l*Nt,(l-1)*Nt+1:l*Nt,l) = sqrt(1/P_BRB(pind))*eye(Nt);
            end
            q = ones(L,1); %Limits of the L power constraints. Note that these have been normalized.
            
            
            %The BRB algorithm searches in a box with lower corner in the
            %origin. The upper corner is the utopia point computed below.
            origin = zeros(Kr,1);
            
            %Computation of the utopia point using MRT, which is optimal
            %when there is only one active user
            utopia = zeros(Kr,1);
            for k = 1:Kr
                utopia(k) = log2(1+abs(H(k,:)*wMRT(:,k))^2*P_BRB(pind));
            end
            
            
            %Obtain the optimal beamforming using the BRB algorithm. Note
            %that this function can produce additional outputs; for example
            %the corresponding beamforming, the number of evaluations of
            %the convex feasibility problem, and the lower and upper bounds
            %when the algorithm terminates.
            bestFeasibleBRB = functionBRBalgorithm_cvx(H,D,Qsqrt,q,origin,utopia,weights,deltaBRB,epsilonBRB,maxIterations,maxFuncEvaluations,bestfeasibleRates(:,pind),problemMode);
            
            %Saving the performance of the optimal beamforming
            sumrateOPTIMAL(pind,m) = weights'*bestFeasibleBRB;
        end
        
    end
    
end




%Plot simulation results
figure; hold on; box on;

if computeOptimalBeamforming==true
    plot(PdB_BRB,mean(sumrateOPTIMAL,2),'k:','LineWidth',1);
end

plot(PdB,mean(sumrateSLNRMAX,2),'r','LineWidth',1);
plot(PdB,mean(sumrateZFBF,2),'b--','LineWidth',1);
plot(PdB,mean(sumrateMRT,2),'k-.','LineWidth',1);

if computeOptimalBeamforming == true
    legend('Optimal','SLNR-MAX','ZF','MRT','Location','NorthWest');
elseif computeOptimalBeamforming == false
    legend('SLNR-MAX','ZF','MRT','Location','Best');
end

xlabel('Average SNR [dB]')
ylabel('Average Sum Rate [bits/symbol use]');
