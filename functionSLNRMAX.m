function wSLNRMAX = functionSLNRMAX(H,eta,D)
%Calculates the Signal-to-leakage-and-noise ratio maximizing (SLNR-MAX)
%beamforming for a scenario where all or a subset of antennas transmit 
%to each user.
%
%This is version 1.0.
%
%INPUT:
%H   = Kr x Kt*Nt matrix with row index for users and column index
%      transmit antennas
%eta = Kr x 1 vector with SNR^(-1) like parameter of this user 
%D   = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if j:th
%      transmit antenna can transmit to user k and zero otherwise
%
%OUTPUT:
%wSLNRMAX = Kt*Nt x Kr matrix with normalized SLNR-MAX beamforming



%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If eta vector is not provided, all values are set to unity
if nargin<2
    eta = ones(Kr,1);
end

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<3
    D = repmat( eye(N), [1 1 Kr]);
end

%Pre-allocation of MRT beamforming
wSLNRMAX = zeros(size(H'));

%Computation of SLNR-MAX, based on Definition 3.5
for k = 1:Kr
    effectivechannel = (H*D(:,:,k))'; %Effective channels
    projectedchannel = (eye(N)/eta(k)+effectivechannel*effectivechannel')\effectivechannel(:,k); %Compute zero-forcing based on channel inversion
    wSLNRMAX(:,k) = projectedchannel/norm(projectedchannel);  %Normalization of zero-forcing direction
end
