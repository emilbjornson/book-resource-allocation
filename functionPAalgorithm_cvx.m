function [bestFeasible,Woptimal,totalNbrOfEvaluations,bounds,vertices] = functionPAalgorithm_cvx(H,D,Qsqrt,q,utopia,weights,delta,epsilon,maxIterations,maxFuncEvaluations,localFeasible,problemMode,saveBoxes)
%Maximizes the weighted sum rate or weighted proportional fairness using
%the Polyblock Outer Approximation (PA) algorithm in Algorithm 2. Both 
%problems are non-convex and NP-hard in general, thus the computational 
%complexity scales exponentially with the number of users, Kr. This 
%implementation is not recommend for more than Kr=6 users.
%
%This is version 1.0.
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%
%INPUT:
%H             = Kr x Kt*Nt matrix with row index for receiver and column
%                index transmit antennas
%D             = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one
%                if j:th antenna can transmit to user k and zero otherwise
%Qsqrt         = N x N x L matrix with matrix-square roots of the L
%              weighting matrices for the L power constraints
%q             = Limits of the L power constraints
%utopia        = Kr x 1 vector with the utopia point, such that the box
%                with lower corner in the origin and upper corner in the
%                utopia point convers the rate region.
%weights       = Kr x 1 vector with positive weights for each user
%delta         = Accuracy of the line-search in FPO subproblems
%                (see functionFairnessProfile() for details
%epsilon       = Accuracy of the final value of the utility
%maxIterations = Maximal number of outer iterations of the algorithm
%maxFuncEvaluations = Maximal number of convex feasibility subproblems to
%                     be solved
%localFeasible = (Optional) Kr x 1 vector with any feasible solution
%problemMode   = (Optional) Weighted sum rate is given by mode==1 (default)
%                 Weighted proportional fairness is given by mode==2
%saveBoxes     = (Optional) Saves and return the set of boxes from each
%                 iteration of the algorithm if saveBoxes==1
%
%OUTPUT:
%bestFeasible          = The best feasible solution found by the algorithm
%Woptimal              = Kt*Nt x Kr matrix with beamforming that achieves bestFeasible
%totalNbrOfEvaluations = Vector with number of times that the convex
%                        subproblem was solved at each iteration of the
%                        algorithm
%bounds                = Matrix where first/second column gives the global
%                        lower/upper bound at each iteration of the algorithm
%vertices              = Cell array where vertices{k}.set contain the
%                        vertices of polyblock at the end of iteration k.



Kr = size(H,1);  %Number of users (in total)


%Initialize the best feasible solution as origin or point given by input
if nargin < 11
    bestFeasible = zeros(Kr,1);
else
    bestFeasible = localFeasible;
end


%If problemMode has not been specified: Select weighted sum rate
if nargin < 12
    problemMode = 1;
end

%If saveBoxes has not been specified: Do not save and return set of boxes
if nargin < 13
    saveBoxes = 0;
end
vertices{1}.set = zeros(Kr,1);


%Initialize origin and vertex set with utopia point
origin = zeros(Kr,1);
vertexSet = max(utopia,[],2);

%Small perturbation term from (2.49) to avoid convergence issues when
%searching for solutions close to an axis.
s = delta/Kr;

%Initialize current best value (cbv) and the current upper bounds (cub),
%where the latter is the potential system utility in each vertex.
if problemMode == 1 %Weighted sum rate
    cbv = weights'*bestFeasible;
    cub = weights'*vertexSet;
elseif problemMode == 2 %Weighted proportional fairness
    cbv = geomean_weighted(bestFeasible,weights);
    cub = geomean_weighted(vertexSet,weights);
end

%Pre-allocation of matrices for storing lower/upper bounds on optimal
%utility and the number of times the convex subproblem (power minimization
%under QoS requirements) is solved.
lowerBound = zeros(maxIterations,1);
upperBound = zeros(maxIterations,1);
totalNbrOfEvaluations = zeros(maxIterations,1);

%Initialize matrix for storing optimal beamforming
Woptimal = zeros(size(H'));


%Iteration of the PA algorithm. Continue until termination by solution
%accuracy, maximum number of iterations or subproblem solutions
for k = 1:maxIterations
    
    %Select vertex with current global upper bound
    [~,refineIndex] = max(cub);
    
    
    %Search on the line between the origin and the selected vertex by
    %solving an FPO problem. This will be used to refine the selected
    %vertex. A feasible point close to the Pareto boundary (accuracy given
    %by delta) is also found and might improve the bounds.
    [interval,W,FPOevaluations] = functionFairnessProfile_cvx(H,D,Qsqrt,q,delta,origin,vertexSet(:,refineIndex));
    
    %Update the number of feasibility evaluations
    totalNbrOfEvaluations(k) = FPOevaluations;
    
    %The new feasible point found by the FPO problem
    newFeasiblePoint = interval(:,1);
    newLocalUpperBound = interval(:,2);
    
    %Update the global feasible point, if the new point is better.
    if problemMode == 1
        if (weights'*newFeasiblePoint) > cbv
            bestFeasible = newFeasiblePoint;
            cbv = weights'*newFeasiblePoint;
            Woptimal = W; %Store beamforming for current best solution
        end
    elseif problemMode == 2
        if geomean_weighted(newFeasiblePoint,weights) > cbv
            bestFeasible = newFeasiblePoint;
            cbv = geomean_weighted(newFeasiblePoint,weights);
            Woptimal = W; %Store beamforming for current best solution
        end
    end
    
    
    
    %Create new vertices based on the new upper bound computed by the FPO
    %problem. The process is based on (2.53) in Lemma 2.8. To avoid
    %numerical convergence issues, new vertices that are very close to an
    %axis are deleted when the distance is smaller than s = delta/Kr
    newVertices = repmat(vertexSet(:,refineIndex),[1 Kr])-diag(vertexSet(:,refineIndex)-newLocalUpperBound);
    notCloseToAnAxis = abs(vertexSet(:,refineIndex))>s;
    newVertices = newVertices(:,notCloseToAnAxis);
    
    if problemMode == 1
        cubNew = weights'*newVertices;
    elseif problemMode == 2
        cubNew = geomean_weighted(newVertices,weights);
    end
    
    
    %Select which vertices to remove. There are three categories: Improper
    %vertices (see last paragraph in Lemma 2.8), Dominated vertices (the
    %local upper bound is lower than the current best feasible solution),
    %and g^(n), which is the vertex that was refined in this iteration.
    removeImproper = (sum(vertexSet>=repmat(newLocalUpperBound,[1 size(vertexSet,2)]),1)==Kr) & (sum(vertexSet<repmat(vertexSet(:,refineIndex),[1 size(vertexSet,2)]))==1);
    removeDominated = cub <= cbv;
    removeOldVertices = removeImproper+removeDominated;
    removeOldVertices(refineIndex) = 1;
    
    
    %Update the vertex set and their local upper bounds for the next 
    %iteration of the algorithm.
    vertexSet = [vertexSet(:,removeOldVertices==0) newVertices];
    cub = [cub(removeOldVertices==0) cubNew];

    
    %Store the lower and upper bounds in the k:th iteration to enable
    %plotting of the progress of the algorithm.
    if problemMode == 1
        lowerBound(k) = cbv;
        upperBound(k) = max(weights'*vertexSet);
    elseif problemMode == 2
        lowerBound(k) = cbv;
        upperBound(k) = max(geomean_weighted(vertexSet,weights));
    end
    
    if saveBoxes == 1
        vertices{k}.set = vertexSet;
    end
    
    
    %Check termination conditions
    if sum(totalNbrOfEvaluations) >= maxFuncEvaluations %Maximal number of feasibility evaluations has been used
        break;
    elseif upperBound(k)-lowerBound(k) <= epsilon %Predefined accuracy of optimal solution has been achieved
        break;
    end
    
end

%Prepare output by removing parts of output vectors that were not used
totalNbrOfEvaluations = totalNbrOfEvaluations(1:k);
bounds = [lowerBound(1:k) upperBound(1:k)];




function y = geomean_weighted(x,w)
%Calculate weighted proportional fairness of each column of z.
%The corresponding weights are given in w

y = prod(x.^repmat(w,[1 size(x,2)]),1);
