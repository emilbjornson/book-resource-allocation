%This m-file generates Figure 3.2 in "Optimal Resource Allocation in
%Coordinated Multi-Cell Systems" by Emil Björnson and Eduard Jorswieck.
%
%This is version 1.2. (Last edited: 2014-07-23)
%
%Figure 3.2 illustrates the performance regions for random channel
%realizations in a scenario with three users and global joint transmission.
%The Pareto boundaries are approximated in two ways: (1) Solving FPO over a
%grid of weighting vectors; and (2) Generating sample points using a
%beamforming parametrization. The color bar shows the sum utility.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the book.

close all;
clear all;


%%Simulation parameters

Kt = 3; %Number of base stations (BSs)
Kr = 3; %Number of users (in total)
Nt = 3; %Number of antennas per BS

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

delta = 0.1; %Accuracy of the line searches in the FPO algorithm


%The grid is created by fixing the performance of one of the users (at
%some different values) and varying the search direction of the FPO for the
%remaining two users.
nbrOfSamplesFPO = 31; %51; %Number of search directions for each fixed value of a user
fixedValuesPerUser = 6;
ParetoBoundary = zeros(Kr,nbrOfSamplesFPO,fixedValuesPerUser,Kr); %Pre-allocation of sample matrix


%Computation of the utopia point using MRT with full transmit power, 
%which is optimal when there is only one active user.

wMRT = functionMRT(H,D); %Compute normalized beamforming vectors for MRT

utopia = zeros(Kr,1);
for k = 1:Kr
    W = sqrt(P)*[wMRT(1:3,k)/norm(wMRT(1:3,k)); wMRT(4:6,k)/norm(wMRT(4:6,k)); wMRT(7:9,k)/norm(wMRT(7:9,k))]; %Scale to use all available power
    utopia(k) = log2(1+abs(H(k,:)*W)^2);
end


%Generate a grid of equally spaced search directions
interval = linspace(0,1,nbrOfSamplesFPO);
profileDirections = [interval; 1-interval];


maxvalue=0.1;

for k = 1:Kr
    
    disp(['Progress: Running Approach 1, User ' num2str(k) ' of ' num2str(Kr)]); %Output of the simulation progress
    
    %Generate the different fixed user rates. The utopia rate is not used,
    %since the other users will achieve very small rates in this case.
    fixedUserRate = linspace(0,utopia(k),fixedValuesPerUser+1); 
    
    for n = 1:fixedValuesPerUser
        lowerPoint = zeros(Kr,1);
        upperPoint = zeros(Kr,1);
        lowerPoint(k) = fixedUserRate(n);
        upperPoint(k) = fixedUserRate(n);
        
        for s = 1:nbrOfSamplesFPO
            
            %Search on a line from the origin in each of the search
            %directions. The upper point is either on Pareto boundary 
            %(unlikely) or outside of the rate region.
            upperPoint([1:k-1 k+1:Kr]) = sum(utopia) * profileDirections(:,s);
            

            %Find the intersection between the line and the Pareto boundary by
            %solving an FPO problem.
            finalInterval = functionFairnessProfile_cvx(H,D,Qsqrt,q,delta,lowerPoint,upperPoint);
            
            %Save the best feasible point found by the algorithm. 
%             if isempty(finalInterval)
%                 ParetoBoundary(:,s,n,k)=lowerPoint;
%             else
            ParetoBoundary(:,s,n,k)=finalInterval(:,1);
%             end
        end
    end
end





%%Approach 2: Calculate sampe points of the whole performance region using
%%the beamforming parameterization in Section 3.2.3.

disp('Progress: Running Approach 2'); %Output of the simulation progress

nbrOfLambdaSamples = 31; %Number of sample points of lambda-parameters
nbrOfMuSamples = 31; %Number of sample points of mu-parameters

lambdaRange=linspace(0,1,nbrOfLambdaSamples);
muRange=linspace(0,1,nbrOfMuSamples);


%Pre-allocated matrices for storing the sample points for the users
feasiblePoints=zeros(nbrOfLambdaSamples^2 * nbrOfMuSamples^2,Kr);

%This counter keeps track of the number of points that have been stored so
%far. This enables storing all points in a 2D matrix.
pointCounter=0;


for lIter1=1:length(lambdaRange)

    for lIter2=1:length(lambdaRange)-lIter1+1
        
        %Sample values on lambda_1, lambda_2, and lambda_3. These have
        %unit sum, as shown by Corollary 3.2.
        
        lambda1 = lambdaRange(lIter1);
        lambda2 = lambdaRange(lIter2);
        lambda3 = 1-lambda1-lambda2;
        
        for mIter1=1:length(muRange)
            for mIter2=1:length(muRange)-mIter1+1
                
                %Sample values on mu_1, mu_2, and mu_3. These have
                %unit sum, as shown by Corollary 3.2.
                
                mu1 = muRange(mIter1);
                mu2 = muRange(mIter2);
                mu3 = 1-mu1-mu2;
                
                %Compute unnormalized versions of the beamforming
                %directions in (3.16). A small number, 1e-8, is added to
                %diagonal to achieve numerical stability.
                w1=(mu1*Q(:,:,1)/q(1)+mu2*Q(:,:,2)/q(2)+mu3*Q(:,:,3)/q(3) + lambda2*H(2,:)'*H(2,:) + lambda3*H(3,:)'*H(3,:) +eye(Kt*Nt)*1e-8)\(H(1,:)');
                w2=(mu1*Q(:,:,1)/q(1)+mu2*Q(:,:,2)/q(2)+mu3*Q(:,:,3)/q(3) + lambda1*H(1,:)'*H(1,:) + lambda3*H(3,:)'*H(3,:) +eye(Kt*Nt)*1e-8)\(H(2,:)');
                w3=(mu1*Q(:,:,1)/q(1)+mu2*Q(:,:,2)/q(2)+mu3*Q(:,:,3)/q(3) + lambda1*H(1,:)'*H(1,:) + lambda2*H(2,:)'*H(2,:) +eye(Kt*Nt)*1e-8)\(H(3,:)');
                
                %Compute gamma_k in (3.19)
                gammas = [real(lambda1*H(1,:)*w1) real(lambda2*H(2,:)*w2) real(lambda3*H(3,:)*w3)];
                
                %Compute normalized beamforming directions as in (3.16)
                V = [w1/norm(w1) w2/norm(w2) w3/norm(w3)];
                
                %Generate M matrix as in (3.20)
                M = zeros(Kr,Kr);
                for i = 1:Kr
                    for j = 1:Kr
                        if i == j
                            M(i,i) = abs(H(i,:)*V(:,i)).^2;
                        else
                            M(i,j) = -gammas(j)*abs(H(j,:)*V(:,i)).^2;
                        end
                    end
                end
                
                %Compute power allocation as in (3.17)
                p = gammas/M;
                
                %Compute final beamforming matrix according to Theorem 3.3
                W = [sqrt(p(1))*w1 sqrt(p(2))*w2 sqrt(p(3))*w3];
                
                %Scale the beamforming matrix according to Corollary 3.2 to ensure feasibility
                maxPowerAmongBSs = max([norm(Qsqrt(:,:,1)*W,'fro').^2 norm(Qsqrt(:,:,2)*W,'fro').^2 norm(Qsqrt(:,:,3)*W,'fro').^2]);
                W = W/sqrt(maxPowerAmongBSs);
                
                
                %Calculate and save user rates with the achieved beamforming
                pointCounter = pointCounter+1;
                channelGains = abs(H*W).^2;
                signalGains = diag(channelGains);
                interferenceGains = sum(channelGains,2)-signalGains;
                feasiblePoints(pointCounter,:) = log2(1+signalGains./(interferenceGains+1))';
                
            end
        end
    end
end



disp('Progress: Final processing for Approach 2'); %Output of the simulation progress

%Go through all sample points and remove all points that are strictly
%dominated by other sample points, since these points cannot be on the 
%Pareto boundary.
for k = 1:pointCounter
    %Check if the m:th point is strictly dominated.
    if sum((feasiblePoints(1:pointCounter,1)>feasiblePoints(k,1)) & (feasiblePoints(1:pointCounter,2)>feasiblePoints(k,2)) & (feasiblePoints(1:pointCounter,3)>feasiblePoints(k,3)))==0
        feasiblePoints(k,:)=0; %Set strictly dominated points to zero 
    end
end

paretoCandidatePoints = feasiblePoints(sum(feasiblePoints,2)>0,:); %Remove strictly dominated points
paretoCandidatePoints = [paretoCandidatePoints; max(paretoCandidatePoints(:,1)) 0 0; 0 max(paretoCandidatePoints(:,2)) 0; 0 0 max(paretoCandidatePoints(:,3)); 0 0 0]; %Add corner points for numerical stability


%Calculate the convex hull of the remaining sample points. This part is the
%main source of computational complexity in Approach 2.
facetsIndices = convhulln(paretoCandidatePoints);





%Plot the rate region generated by Approach 1.

figure(1);

plot3(ParetoBoundary(2,:,1,1),ParetoBoundary(1,:,1,1),ParetoBoundary(3,:,1,1),'k'); %Plot first curve to get the right dimensions

hold on; box off;

%Iteration over all grid points
for k = 1:Kr
    for n = 1:fixedValuesPerUser
        plot3(ParetoBoundary(2,:,n,k),ParetoBoundary(1,:,n,k),ParetoBoundary(3,:,n,k),'k');
    end
end
view([120 30]);
axesApproach1 = axis;

ylabel('log_2(1+SINR_1)')
xlabel('log_2(1+SINR_2)')
zlabel('log_2(1+SINR_3)')


%Plot the rate region generated by Approach 2.

colors=sum(paretoCandidatePoints,2); %The color of each sample point corresponds to the sum rate

figure(2); grid off; box off;

trisurf(sort(facetsIndices,2),paretoCandidatePoints(:,2),paretoCandidatePoints(:,1),paretoCandidatePoints(:,3),colors,'EdgeColor','none','FaceColor','interp','CDataMapping','scaled')
view([120 30]);
colormap(flipud(hot)); colorbar; %The color bar shows the sum rate

ylabel('log_2(1+SINR_1)')
xlabel('log_2(1+SINR_2)')
zlabel('log_2(1+SINR_3)')

axis(axesApproach1);
