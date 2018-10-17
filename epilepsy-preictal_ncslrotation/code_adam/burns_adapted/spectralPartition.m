%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script spectralPartition.m
%
% Description: Use the spectral partitioning algorithm to cluster the
% EVC's by making a network of the EVC's in interictal, preictal, ictal and
% postictal over time and undergoing the spectral analysis.
%
% Input: the graph produce by 'evcNetwork.m' 
%
% Output: the different clusters that each EVC in time belongs to
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 08/21/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [cluster1 cluster2] = spectralPartition(adjMat)
    %% 01: Create Settings for Running Algorithm
    homedir = pwd;
    
    adjMat = graph;
    % compute diagonal vertex matrix
    degreeMat = diag(sum([adjMat adjMat'],2)) - diag(diag(adjMat));
    
    % compute Laplacian(adjMat)
    L = degreeMat - adjMat;
    [V, D] = eigs(L);
    
    %% 02: compute Fiedler Eigenvector
    FiedlerEV = V(:,2); % get the eigenvector for second largest eigenvalue
    sortedEV = sort(FiedlerEV);
    lambda2 = D(2);     % store second largest eigenvector
    
    %% 03: Compute median(Fiedler Vector)
    medV = median(FiedlerEV)
    
    vertex1 = [];
    vertex2 = [];
    %% 04: For each node 'i' in adjacency matrix, put node in either V1, or V2 
    for i=1:length(FiedlerEV)
       if FiedlerEV(i) <= medV
           vertex1(end+1) = i;
       elseif FiedlerEV(i) > medV
           vertex2(end+1) = i;
       end
    end


    %% 05: Check if the count of  cluster1 == cluster2
    if abs(length(vertex1) - length(vertex2)) > 1
       disp('Many vertices equal to median! Fix Step 05:')
    end
    
    % calculate cut size R = n1*n2/n * lambda
    Rcutsize = (length(vertex1) * length(vertex2)) / (length(FiedlerEV)) * lambda2
    %% 06: Get V1' and V2' the set of edges in sets V1, V2 adj. to some vertex in V1, V2
    edge1 = [];
    edge2 = [];
    for i=1:length(vertex1)
        for j=1:length(vertex1)
            
        end
    end
    
    for i=1:length(vertex1)         % sum over vertex1 group
        for j=1:length(vertex2)     % sum over vertex2 group
            
        end
    end
    
% end