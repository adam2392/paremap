function [allmembersCluster clustStatsOut clustStatsInd] = getSigClustB_clustStats(data,adj,sigThresh,clustStatList,clustLimit,minClustSize)
% [sumCluster,membersCluster,allmembersCluster] = getSigClust_clustStats(data,adj)
%
%  INPUT:
%    data (numROI x numPerms)     = thresholded but weighted effect matrix (0 or analog value of "significant" ROI.  can be pos or neg)
%    adj  (numROI x numROI)       = binary symmetric unweighted adjacency matrix of ROIs
%    sigThresh                    = value used to threshold the data.  data should already be thresholded -- this is only used for computing cluster Mass
%    clustStatList {}             = select cluster stat(s) with cell of strings:  'size' 'sum' 'max' 'mass'.  empty cell array means dont compute cluster stats (still find allmembers and return that)
%    clustLimit                   = enter number of TRUE clusters here to constrain clustStatOut dimension. Set to -1 to include all clusters from each perm (use that for "TRUE" case)
%    minClustSize (>=1)           = minimum cluster size for considering cluster "real" and calculating cluster statistics.  for ROI's use 1, for Spectrograms use ~5 or 10. 
%
%  OUTPUTS:
%    allmembersCluster {numPerms} = cell of cells: each perm's cell --> cell array of cluster members
%    clustStatsOut {numStats}([numPerm x numclust]) = cell of matricies: each clustStat gets a cell entry, that entry is a numPerm x clustLimit matrix of sorted cluster stats
%    clustStatsInd {numStats}([numclust])           = cell of arrays: index of allmemberClusters that map to sorted StatsOut. Only the last perm's StatsInd is saved (only need this for true, which has 1 perm)
%
%  JHW 1/2016 (derived from stuff Rafi and Vaz gave me)
%

%- this should be 1 for ROI, 4? for spectrograms... or only apply the trim to the "TRUE" clusters in getSigClustA? seems more fair to apply to both here
MIN_CLUS_SIZE = minClustSize; 
if MIN_CLUS_SIZE>1 & clustLimit>0, fprintf('\n HEADS UP: eliminating all clusters < %d from true and/or perm \n',MIN_CLUS_SIZE); end


% Reshape matrix (e.g. numROIs x numPerms), you will cluster each t,f,p combination
dataSize = size(data);  
numPerms = prod(dataSize(2:end));
if length(dataSize)>2, fprintf('\n Uh oh.. only setup for ROI x nPerm but >2 dimensions passed in'); keyboard; end


% Initialize cluster stat matricies and identify the (potential) output order so that it matches "clustStatList"
getSize = 0; getSum  = 0; getMax  = 0; getMass = 0; getAny = 0;
NUM_CLUST_OUT = clustLimit; if NUM_CLUST_OUT<0, NUM_CLUST_OUT=0; end;  %- use this for initializing perm matricies
iSize = find(strcmp(clustStatList,'size'));  if ~isempty(iSize), getSize = 1;  perm_clusSiz = zeros(numPerms,NUM_CLUST_OUT); iSrtSiz = []; end %- size of each cluster
iSum  = find(strcmp(clustStatList,'sum'));   if ~isempty(iSum),  getSum  = 1;  perm_clusSum = zeros(numPerms,NUM_CLUST_OUT); iSrtSum = []; end %- sum  = sum of z-scores of all nodes in cluster
iMax  = find(strcmp(clustStatList,'max'));   if ~isempty(iMax),  getMax  = 1;  perm_clusMax = zeros(numPerms,NUM_CLUST_OUT); iSrtMax = []; end %- max  = maximum value of all nodes in cluster 
iMass = find(strcmp(clustStatList,'mass'));  if ~isempty(iMass), getMass = 1;  perm_clusMas = zeros(numPerms,NUM_CLUST_OUT); iSrtMas = []; end %- mass = sum of values above threshold (could be zero if all nodes are exactly at threshold)
if NUM_CLUST_OUT==0, NUM_CLUST_OUT=inf; end;  %- this means all cluster info will be saved
getAny = getSize+getSum+getMax+getMass;
clustStatsOut = cell([1 getAny]);
clustStatsInd = cell([1 getAny]);

allmembers    = cell([1 numPerms]);  %- always get all the members, even if no clust stat is computed


% Loop through each permutation and find the clusters + cluster statistics
tic;
if numPerms>1, divPerms=round(numPerms/10); else divPerms=1; end
for iPerm = 1:numPerms,
    % output progress to commandline 
    if mod(iPerm,divPerms)==1, fprintf('%.0f%% ',round(100*iPerm/numPerms)); end  
    
    % access specific perm you want to cluster
    sigData = data(:,iPerm)';
    iSig    = find(sigData);  %- indices to non-zero entries, these are the only ones to bother trying to cluster 
      
    if length(iSig)>0, 
        %-- at least 1 significant ROI --%
        
        % remove ROIs that aren't significant (doesn't matter if their neighbors are)
        newClusterThis = adj(iSig,iSig);  
        
        % find clusters
        [tmpMembers] = getSigClustC_findClust(newClusterThis);  %- clusters based on unweighted matrix (treats as 0 and nonzero). every node belongs to a cluster(possibly self), even if that node isn't sig... so only pass sig nodes
       
        % cut small clusters here... means they wont even show up in "uncorrected" node output
        if MIN_CLUS_SIZE>1,  tmpMembers=tmpMembers(find(cellfun(@length, tmpMembers)>=MIN_CLUS_SIZE)); end
            
        % save cluster members (first convert back to original ROI space).  sig cluster can contain as few as 1 ROI that is significant but not connected
        allmembers{iPerm} = cellfun(@(x) iSig(x),tmpMembers,'uniformoutput',false); %-cell array: list of ROIs (translated to original ROI space using iSig)

        
        %- compute the cluster stats, sort, and save up to NUM_CLUST_OUT
        if getAny,
            
            %- For spectrograms only consider clusters >= MIN_CLUS_SIZE for multivariate cluster correction
            %             if MIN_CLUS_SIZE>1,
            %                 iTmpMem=find(cellfun(@length, tmpMembers)>=MIN_CLUS_SIZE); tmpMembers=tmpMembers(iTmpMem);
            %             else
            %                 iTmpMem=1:length(tmpMembers); %-
            %             end
            iTmpMem=1:length(tmpMembers);
            
            iClus = 1:min([length(tmpMembers) NUM_CLUST_OUT]);
            
            
            %- Note: following iSrt indicies only important for "true" perm so contributing nodes can be identified, dont bother saving in perm x clus matrix (just clust array)
            %-  size of cluster  
            if getSize, [sVals,iSrtSiz] = sort(cellfun(@length,                                   tmpMembers),'descend'); iSrtSiz=iTmpMem(iSrtSiz); perm_clusSiz(iPerm,iClus) = sVals(iClus);  end 
            %-  absolute value of sum(Z) [pos and neg in same cluster will cancel]
            if getSum,  [sVals,iSrtSum] = sort(cellfun(@(x) abs(sum(sigData(iSig(x)))),           tmpMembers),'descend'); iSrtSum=iTmpMem(iSrtSum); perm_clusSum(iPerm,iClus) = sVals(iClus);  end 
            %-  max absolute Z value
            if getMax,  [sVals,iSrtMax] = sort(cellfun(@(x) max(abs(sigData(iSig(x)))),           tmpMembers),'descend'); iSrtMax=iTmpMem(iSrtMax); perm_clusMax(iPerm,iClus) = sVals(iClus);  end 
            %-  sum absolute excees Zthresh [pos and neg is same cluster will NOT cancel]
            if getMass, [sVals,iSrtMas] = sort(cellfun(@(x) sum(abs(sigData(iSig(x)))-sigThresh), tmpMembers),'descend'); iSrtMas=iTmpMem(iSrtMas); perm_clusMas(iPerm,iClus) = sVals(iClus);  end 
        end
       
    else
        %-- no significant Nodes or Clusters found --%
        allmembers{iPerm} = [];
    end
end
if numPerms>1, fprintf(' [%.1f s]',toc); end

% reshape outputs into the sizes of inputs
numIndFeatSize    = dataSize(2:end);
allmembersCluster = reshape(allmembers,[1 numIndFeatSize]);

if getSize, clustStatsOut{iSize} = perm_clusSiz; clustStatsInd{iSize} = iSrtSiz; end
if getSum,  clustStatsOut{iSum}  = perm_clusSum; clustStatsInd{iSum}  = iSrtSum; end
if getMax,  clustStatsOut{iMax}  = perm_clusMax; clustStatsInd{iMax}  = iSrtMax; end
if getMass, clustStatsOut{iMass} = perm_clusMas; clustStatsInd{iMass} = iSrtMas; end

