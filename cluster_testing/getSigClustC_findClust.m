% [members] = getSigClustC_findClust(A)
%
% Daniel Larremore
% May 31, 2013
% larremor@hsph.harvard.edu
% http://danlarremore.com
% Comments and suggestions always welcome.
%
% INPUTS:
% A                     Matrix. fThis function takes as an input a
% network adjacency matrix A, for a network that is undirected. If you
% provide a network that is directed, this code is going to make it
% undirected before continuing. Since link weights will not affect
% component sizes, weighted and unweighted networks work equally well. You
% may provide a "full" or a "sparse" matrix.
%
% OUTPUTS:
% members               cell<vector<INT>> a cell array of vectors, each
%   entry of which is a membership list for that component.%
%
%  1/2016  JW version 2: only return members to make a little faster... cluster stats computed next level up
%          could convert this to expect a sparce cell array adjacency matrix: just a list of nodes connected to N... ah, but this would make a trimmed version tricky, would have to nan non-sig values

function [members] = getSigClustC_findClust(A)

% Number of nodes
N = size(A,1);
% Remove diagonals
A(1:N+1:end) = 0;
% make symmetric, just in case it isn't
A=A+A';
% Have we visited a particular node yet?
isDiscovered = zeros(N,1);
% Empty members cell
members = {};

% check every node
for n=1:N
    if ~isDiscovered(n)
        % started a new group so add it to members
        members{end+1} = n;
        % account for discovering n
        isDiscovered(n) = 1;
        % set the ptr to 1
        ptr = 1;
        while (ptr <= length(members{end}))
            % find neighbors
            nbrs = find(A(:,members{end}(ptr)));
            % here are the neighbors that are undiscovered
            newNbrs = nbrs(isDiscovered(nbrs)==0);
            % we can now mark them as discovered
            isDiscovered(newNbrs) = 1;
            % add them to member list
            members{end}(end+1:end+length(newNbrs)) = newNbrs;
            % increment ptr so we check the next member of this cluster
            ptr = ptr+1;
        end
    end
end

% % number of clusters
% nClusters = length(members);
% for n=1:nClusters
%     % compute sizes of clusters
%     sizes(n) = length(members{n});
% end
% % find the biggest one
% idxLargestCluster = find(sizes==max(sizes));
% end