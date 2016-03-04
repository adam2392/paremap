function [allClust_uncor, sigClust_stat, sigClust_ind, sigClust_pVal, strOut] = getSigClustA_mainROIs(truedata,permdata,adj,sigThresh,lumpPosNeg,minClustSize)
% function getSigClustA_mainROIs -- execute cluster correction.  multiple options for cluster statistic and analysis
%
% INPUTS:
%   truedata (ROI x 1)
%   permdata (ROI x numPerm)
%   adj      (ROI x ROI)  - adjacency matric
%   sigThresh  (double)   - threshold value used to convert true and perm data into 0 and sig values
%   lumpPosNeg (0 or 1)   - (0)=adjacent sig ROIs of same sign cluster; (1)=adjacent significant ROIs cluster irrespective of sign;  
%   minClustSize (>=1)    - minimum cluster size for considering a cluseter "real" and calculating cluster statistics.  for ROI's use 1, for Spectrograms use ~5 or 10
%
% OUTPUTS:
%   allClust_uncor {numClust}(numTFxthisClust) = list of inidcies into truedata that make up each cluster. if numTF>1, it will be a 2-row output, with ROI on the top and TimeFreq index on the bottom
%   sigClust_stat  {numStats}              = cell array of strings defining each cluster stat and correction method used (e.g.,  'size: vs biggest', 'size: stepDown', 'mass: multiVar FDR', 'mass: multivar Liptak')
%   sigClust_ind   {numStats}(numSigClust) = list of indicies of significant clusters (references allClust list) 
%   sigClust_pVal  {numStats}(numSigClust) = list of p-values of significant clusters (references allClust list) 
%
%
%  JHW 1/2016 (derived from stuff Rafi and Vaz gave me)
%


%%------- SET ANALYSIS TOGGLES HERE -------%
whichStat = {'sum','size','max','mass'};   %- {'size','sum','max','mass'}:  SELECT 0,1, or more cluster stats
USE_STEP_DOWN_METHOD     = 0;  %%- STEP-DOWN METHOD:  if biggest true clust is significant, remove from true and permuted data and recompute perm clusters to find next biggest (iterate until not sig)
USE_MATCHED_CLUST_METHOD = 1;  %%- MATCHED CLUST METHOD: MULTIVARIATE analysis of clusters, where each cluster is compared to matched-rank cluster in permuted distribution
matchedClustTypes = {'FDR','Fisher','Liptak','Tippett'}; %- matchedClustTypes = {'FDR','Fisher','Liptak','Tippett'};    SELECT at least 1 matchedClusterType if USE_MATCHED_CLUST_METHOD=1.  

MAX_CLUST_LIPTAK         = 15; %%- this is a computational limitation:  office computer can do 20 in 30 seconds; 18 in 7 seconds; 15 in .6 seconds
        
whichStat         = {'sum'};   %- {'size','sum','max','mass'}:  SELECT 0,1, or more cluster stats
matchedClustTypes = {'FDR'};   %- matchedClustTypes = {'FDR','Fisher','Liptak','Tippett'};    SELECT at least 1 matchedClusterType if USE_MATCHED_CLUST_METHOD=1.  




%- initialize variables and output text -%
clustStatFig = 0;
sigClust_stat={}; sigClust_ind={[]}; sigClust_pVal={[]}; strOut = '';
numROI = size(permdata,1); numPerm=size(permdata,2); numTF=size(permdata,3);
fprintf('\n Cluster Correction (%d ROI x %d PERMS x %d Time-OR-Freqs):',numROI,numPerm,numTF); tStart=tic;


%- if more than one time OR frequency point reshape the data into a long list of data values and expand the adjacency matrix so that sequential ROI x perm matricies are considered adjacent
%    NOT valid if time AND frequency are multiplexed.. that would require an outer loop to the iTF loop below
if numTF>1,
    truedata = reshape(permute(truedata,[1 2]),numROI*numTF,1);            %- stack true and perm time or freq data so its one long vector/matrix
    permdata = reshape(permute(permdata,[1 3 2]),numROI*numTF,numPerm);
    adj2 = zeros(numROI*numTF);
    for iTF=1:numTF,
       iAdj = [1:numROI]+numROI*(iTF-1) ;                                  %- replicate the spatial map for each time/frequency point
       adj2(iAdj,iAdj) = adj;
       if iTF>1,     adj2(iAdj,iAdj-numROI) = diag(ones(numROI,1)); end    %- each time or freq point is adjacent to itself in the neighboring spatial map
       if iTF<numTF, adj2(iAdj,iAdj+numROI) = diag(ones(numROI,1)); end
    end
    %figure(1); clf;  %- sanity check that folding out the time/freq points is correct
    %subplot(221); imagesc(adj); set(gca,'ydir','normal'); axis equal; axis square; 
    %subplot(223); imagesc(diag(ones(numROI,1))); set(gca,'ydir','normal'); axis equal; axis square;
    %subplot(122); imagesc(adj2); set(gca,'ydir','normal'); axis equal; axis square;    
    adj = adj2;
end


%- threshold and perform clustering on the true case... this determines the size of the cluster-statistic vector collected from permutations
truedata(isnan(truedata)) = 0;
if lumpPosNeg==1,
    strLump = 'Lump';
    truedata(abs(truedata) < sigThresh) = 0; %-pos and neg clustered together
    [allmembers_true,clustStats_true,clustStatInd_true] = getSigClustB_clustStats(truedata,adj,sigThresh,whichStat,-1,minClustSize);  %allmembersCluster clustStatsOut
else
    strLump = 'Sep';
    trueBoth = truedata;
    truedata(  truedata  < sigThresh) = 0;   %-pos
    [allmembers_truePos,clustStat_truePos,clustStatInd_truePos] = getSigClustB_clustStats(truedata,adj,sigThresh,whichStat,-1,minClustSize);
    truedata = trueBoth;
    truedata( -truedata  < sigThresh) = 0;   %-neg
    [allmembers_trueNeg,clustStat_trueNeg,clustStatInd_trueNeg] = getSigClustB_clustStats(truedata,adj,sigThresh,whichStat,-1,minClustSize);
    %-recombine data so it can all be evaluated at once
    truedata = trueBoth;
    allmembers_true{1} = [allmembers_truePos{:} allmembers_trueNeg{:}];
    for iStat=1:length(whichStat),
        %- re-sort the merged positive and negative sets
        [clustStats_true{iStat},sortInd] = sort([clustStat_truePos{iStat}    clustStat_trueNeg{iStat}],'descend');
        %- update the index array into allmembers so the sorted clusters are properly associated with ROI locations
        clustStatsInd_true{iStat}        = [clustStatInd_truePos{iStat} clustStatInd_trueNeg{iStat}+length(allmembers_truePos{:})];
        clustStatsInd_true{iStat}        = clustStatsInd_true{iStat}(sortInd);
    end
end

NUM_TRUE_CLUST = length(allmembers_true{1}); %- already trimmed based on size for spectrograms
allClust_uncor = allmembers_true{1}; %- a row vector of ROI indicies for each cluster;  if numTF>1, then it will become a 2-row matrix of ROI and TimeFreq indicies
if numTF>1,  
    for iClust = 1:NUM_TRUE_CLUST,
        nodeList = allClust_uncor{iClust};         %- nodeList = [1:numROI*numTF]
        roiList  = mod((nodeList-1),numROI)+1;   
        tfList   = floor((nodeList-1)/numROI)+1;
        %fprintf('\n %04d  \t  %03d  \t  %02d',[nodeList' roiList' tfList']'); %- linearized ROIxTF, ROI, TF
        
        allClust_uncor{iClust} = [roiList' tfList']'; % create 2D outputs that point to ROI index and Time-OR-Freq index
    end
end


%- catch here... if no true clusters, don't bother with permutation
if NUM_TRUE_CLUST==0,
    fprintf(' no true clusters, so jump out of this function');  
    sigClust_stat={'no true clust'};
    if numPerm==0,
        figure(81); clf; %- so plot from below is erased
    end
    return;
    
elseif numPerm==0,
    
    %- figure of clust stats
    figure(81); clf;  set(gcf,'position',[350  800  2000  350]);
    for iStat=1:length(whichStat),
        trueStat = clustStats_true{iStat};
    
        %- cluster stats
        subplot(1,length(whichStat),iStat);
        plot([1:NUM_TRUE_CLUST],  trueStat, 'r*','markersize',20);  hold on; grid on;
        set(gca,'fontsize',15,'ylim',[0 max(trueStat)*1.1],'xlim',[0 NUM_TRUE_CLUST+1],'box','off');
        title(sprintf('Cluster %s',whichStat{iStat}));
        xlabel('Cluster Number');
        if iStat==1, ylabel('True Cluster Stat'); end
        if iStat==1, text(1,max(trueStat)*1.05,sprintf('[ROIxTFxPERMS= %d x %d x %d; min clust size=%d]',numROI,numTF,numPerm,minClustSize),'fontsize',16); end
        if strcmp(whichStat{iStat},'size'), set(gca,'yscale','log'); end
    end
    drawnow;
    return;
    
elseif USE_MATCHED_CLUST_METHOD,
    NUM_PERM_CLUST_SAVE = NUM_TRUE_CLUST;
else
    NUM_PERM_CLUST_SAVE = 1;
end



%- threshold and perform clustering on the permuted cases
permdata(isnan(permdata)) = 0;
if lumpPosNeg==1, permdata(abs(permdata) < sigThresh) = 0;     %-pos and neg clustered together... this can create bigger cluster that have smaller sum values
else              permdata(    permdata  < sigThresh) = 0; end %-pos and neg separate -- for across subject perms emperic distribution is symmetric, so only need to compute 1 (if compute then combine dist is identical)
[allmembers_perm,clustStatsOut_perm,~] = getSigClustB_clustStats(permdata,adj,sigThresh,whichStat,NUM_PERM_CLUST_SAVE,minClustSize);
% INPUTS:  permdata should be thresholded and weighted (0 or significant non-zero value); adj is ROI-2-ROI adjacency matrix
% OUTPUTS: allmembers_perm{nPerm}{allClustPerPerm}; clustStatsOut_perm{numStats}(numPerm,NUM_TRUE_CLUST)


%- tweak these if no perms...??
%if numPerm==0, USE_STEP_DOWN_METHOD=0; USE_MATCHED_CLUST_METHOD=0;  end

%- Loop through the cluster statistics and evaulate for significance (size, sum, max, or mass)
numStats = length(whichStat);
for iStat=1:numStats,
    
    clustStr = sprintf('Cluster %s',whichStat{iStat});
    permStat = clustStatsOut_perm{iStat};
    trueStat = clustStats_true{iStat};
    
    %- get a pValue for each cluster vs the biggest cluster stat in each permutation (single test controls for family wise error)
    for iClust = 1:NUM_TRUE_CLUST,
        truePvsBiggest(iClust) = (sum( permStat(:,1) >= trueStat(iClust) )+0.5) / (numPerm+1);  %- one tailed:  p = t *>= t  (only looking for biggest, not biggest and smallest)
    end
    sigClust_stat{iStat} = sprintf('%s: vs biggest',whichStat{iStat});
    sigClust_ind{iStat}  = clustStatsInd_true{iStat}(truePvsBiggest<=0.05);
    sigClust_pVal{iStat} = truePvsBiggest(truePvsBiggest<=0.05);
    
    
    %%- STEP-DOWN METHOD:  if biggest true clust is significant, remove from true and permuted data and recompute perm clusters to find next biggest (iterate until not sig)
    if USE_STEP_DOWN_METHOD,
        %- initialize pStepDown = truePvsBiggest... then fill in subsequent values with each perm
        truePstepDown = truePvsBiggest;
        iClust  = 1;
        roiKeep = [1:numROI];
        
        %- while biggest true is significant, remove from dataset and repermute to test 2nd biggest true
        while truePstepDown(iClust)<=0.05,
           trueClustInd = allmembers_true{1}{clustStatsInd_true{iStat}(iClust)}; %- ROIs in true significant cluster
           roiKeep = setdiff(roiKeep,trueClustInd);
           fprintf('\n clust %d sig, cutting %d ROI and re-permuting (%d of %d ROI left): ',iClust,length(trueClustInd),length(roiKeep),numROI);
            
           [~,thisPermStat,~] = getSigClustB_clustStats(permdata(roiKeep,:),adj(roiKeep,roiKeep),sigThresh,whichStat{iStat},1); %- just save max for stepdown
           iClust=iClust+1;
           thisPermStat = thisPermStat{1};
           truePstepDown(iClust) = (sum( thisPermStat(:,1) >= trueStat(iClust) )+0.5) / (numPerm+1);  %- one tailed:  p = t *>= t  (only looking for biggest, not biggest and smallest)
        end
        
        sigClust_stat{iStat} = sprintf('%s: step down',whichStat{iStat});
        sigClust_ind{iStat}  = clustStatsInd_true{iStat}(truePstepDown<=0.05);
        sigClust_pVal{iStat} = truePstepDown(truePstepDown<=0.05);
    end
    
    
    %%- MULTIVARIATE analysis of clusters, where each cluster is compared to matched-rank cluster in permuted distribution
    if USE_MATCHED_CLUST_METHOD,
        
        fprintf(' -%s- ',clustStr);
    
        %- get a pValue for each cluster vs the matched-rank cluster from each permutation (biggest vs biggest, 2nd vs 2nd, etc)
        %    requires followup correction for family wise error because mulitple comparisons
        combPval = mat2pval([permStat; trueStat]);  %- calculate the permutation and true p values at the same time
        permPval = combPval(1:end-1,:);
        truePvsMatched = combPval(end,:); %- this is legit... possible to be <0.05 if true > 10 of 10 perms (0/10 >= true --> p=0.5/11)
        
        
        %- Use NPC package to effeciently compute non-parametric combining function on permuted data, corrected for Family Wise Error;
        % Fisher  combining function:  sum of -log(p) --> intermediate
        % Liptak  combining function:  sum of Z's     --> best for weak global
        % Tippett combining function:  max p          --> best for localized    fprintf('   %s: ',clustStr);
        %   combPval = [permPval; truePvsMatched];  %--> NPC software expects last entry to be true pvalues
        if NUM_TRUE_CLUST>MAX_CLUST_LIPTAK & (sum(strcmp(matchedClustTypes,'Fisher')) | sum(strcmp(matchedClustTypes,'Liptak'))),
            if iStat==1, fprintf('\n %d true clusters; only evaluating first %d for Fisher and Liptak: ',NUM_TRUE_CLUST,MAX_CLUST_LIPTAK); end
            adj_P_Fisher = nan*truePvsMatched;
            adj_P_Liptak = nan*truePvsMatched;
            cRng = 1:MAX_CLUST_LIPTAK;
        else
            cRng = 1:NUM_TRUE_CLUST;
        end
        
        if sum(strcmp(matchedClustTypes,'Fisher')),
            tic; [adj_P_Fisher(cRng), p_glob_Fisher, ~] = NPC_FWE_jw(combPval(:,cRng),'F'); fprintf('F %.1fs; ',toc);  % Fisher  combining function:  sum of -log(p) --> intermediate
            adj_P_out = adj_P_Fisher;  strAdjType = 'Fisher';
        end
        if sum(strcmp(matchedClustTypes,'Liptak')),
            tic; [adj_P_Liptak(cRng), p_glob_Liptak, ~] = NPC_FWE_jw(combPval(:,cRng),'L'); fprintf('L %.1fs; ',toc);  % Liptak  combining function:  sum of Z's     --> best for weak global
            adj_P_out = adj_P_Liptak;  strAdjType = 'Liptak';
        end
        if sum(strcmp(matchedClustTypes,'Tippett')),
            tic; [adj_P_Tippet, p_glob_Tippet, ~] = NPC_FWE_jw(combPval,'T'); fprintf('T %.1fs; ',toc);                % Tippett combining function:  max p          --> best for localized
            adj_P_out = adj_P_Tippet;  strAdjType = 'Tippett';
        end
        
        
        %- Cheap and quick FDR of multiple cluster p-values
        if sum(strcmp(matchedClustTypes,'FDR')),
            [pThreshFDR,~] = fdr(truePvsMatched,0.05); if isempty(pThreshFDR), pThreshFDR=nan; end
            [adj_P_FDR]    = pAdjust(truePvsMatched);
            adj_P_out = adj_P_FDR;  strAdjType = 'FDR';
        end
        
        
        %- looks like Liptak is the most inclusive for brain plots... Tippett for spectrograms where LOTS of clusters exist. FDR might work for everything
        sigClust_stat{iStat} = sprintf('%s: %s',whichStat{iStat},strAdjType);
        sigClust_ind{iStat}  = clustStatsInd_true{iStat}(adj_P_out<=0.05);
        sigClust_pVal{iStat} = adj_P_out(adj_P_out<=0.05);
    end
    
    
    %- plot the ranked cluster statistics (true vs uncorrected 95-th percentile)
    if (1 & length(matchedClustTypes)==4),
        
        
        figure(80); if iStat==1, clf; thisAx=[]; set(gcf,'position',[700    200   1800  900]); drawnow; pause(1); end; %- drawnow here to eliminate problem where column one headings passed back to calling function?
        iX = [1:NUM_TRUE_CLUST];
        
        %- cluster stats
        subplot(3,numStats,iStat); thisAx=[thisAx gca];
        plot(iX, trueStat, 'r*','markersize',20);  hold on; 
        set(gca,'fontsize',15,'box','off'); title(clustStr);
        plot([0 NUM_TRUE_CLUST+1], [1 1]*quantile(permStat(:,1),.95), 'k:');
        plot(iX, quantile(permStat,.95), 'k--');
        xlabel('Cluster Number');
        if iStat==1, ylabel('Statistic');  legend('true cluster stat','95th percentile biggest','uncorrected 95th percentile','location','best'); end
        
        
        %- p-values
        subplot(3,numStats,iStat+numStats); thisAx=[thisAx gca];
        semilogy([0 NUM_TRUE_CLUST+1],[1 1]*0.05, 'k:','linewidth',3); hold on; 
        set(gca,'fontsize',15, 'box','off'); grid on;
        h0 = semilogy(iX, truePvsBiggest,'r.','markersize',15);
        if USE_STEP_DOWN_METHOD,
            hX = semilogy(iX, truePstepDown,'rs','markersize',10);
        end
        if USE_MATCHED_CLUST_METHOD,
            semilogy([0 NUM_TRUE_CLUST+1],[1 1]*pThreshFDR,'k:','linewidth',2);
            h1 = semilogy(iX, truePvsMatched,'r*','markersize',15);
            h2 = plot(iX+.1,  adj_P_Fisher,'g*','markersize',6);
            h3 = plot(iX-.1,  adj_P_Liptak,'m*','markersize',6);
            h4 = plot(iX+.05, adj_P_Tippet,'b*','markersize',6);
            if iStat==1, legend([h0 h1 h2 h3 h4],'p vs biggest','Uncorrected p','FWE Fisher','FWE Liptak','FWE Tippet','location','best'); end
        end
        if iStat==1, ylabel('uncorrected pvalue'); end
        
        
        %- significant or not
        subplot(3,numStats,iStat+numStats*2); thisAx=[thisAx gca]; 
        iSig = truePvsBiggest<=0.05; ySig=4; h0=plot([-1 iX(iSig)], ySig, 'r*', 'markersize',15);  hold on;   %- plot dummy point at -1 so legend can always be populated
        linkaxes(thisAx,'x');  %- link first so following xlim change affects everything
        set(gca,'fontsize',15,'ylim',[-0.5 4.5],'yticklabel',{},'xlim',[0 NUM_TRUE_CLUST+1],'box','off')
        if USE_STEP_DOWN_METHOD,
            iSig = truePstepDown<=0.05; ySig=4; hX=plot([-1 iX(iSig)], ySig, 'rs', 'markersize',15);  hold on;   %- plot dummy point at -1 so legend can always be populated
        end
        if USE_MATCHED_CLUST_METHOD,
            iSig = adj_P_Fisher<=0.05; ySig=3; h1=plot([-1 iX(iSig)], ySig, 'g*', 'markersize',15);  hold on;   %- plot dummy point at -1 so legend can always be populated
            iSig = adj_P_Liptak<=0.05; ySig=2; h2=plot([-1 iX(iSig)], ySig, 'm*', 'markersize',15);
            iSig = adj_P_Tippet<=0.05; ySig=1; h3=plot([-1 iX(iSig)], ySig, 'b*', 'markersize',15);
            iSig = adj_P_FDR   <=0.05; ySig=0; h4=plot([-1 iX(iSig)], ySig, 'k*', 'markersize',15);
            if iStat==1, legend([h0(1) h1(1) h2(1) h3(1) h4(1)],'Biggest','Fisher','Liptak','Tippett','FDR','location','best'); end
        end
        
        if iStat==1, ylabel('Significant Cluster (corrected');  xlabel(sprintf('[%d perms; %d TF; min clust size (true and perm)=%d]',numPerm,numTF,minClustSize)); end
        drawnow;
    end
end




PLOT_RAW_CLUSTER_DIST=0;
if PLOT_RAW_CLUSTER_DIST,
    
    %- plot the raw cluster stat distributions ---
    
    fprintf('\n following code isnt right yet... needs to finish being updated to "clustStatsOut_perm" and "whichStat" notation');
    keboard
    
    figure(69); clf;
    clusCount = max([5 NUM_TRUE_CLUST]);  axList = [];
    
    numStats = length(whichStat);
    for iStat=1:numStats,
        
        for iClus=1:clusCount,
            axList(iClus,1) = subplot(clusCount,numStats,1+(iClus-1)*numStats);
            hist(perm_clusSiz(:,iClus),[0:50]); axis tight; hold on; box off; set(gca,'fontsize',15)
            plot([1 1]*quantile(perm_clusSiz(:,iClus),.95),get(gca,'ylim'),'r--','linewidth',2); hold on;
            if iClus==1, title('Size'); legend('emperic dist','95th pcntile'); end
            
        end
        linkaxes(axList(1:clusCount,iStat),'x');
    end
end



%fprint('\n OUTPUTS NEED TO BE CLEANED UP');


%   sigUncor_nodes (num sig nodes)         = list of inidcies into truedata
%   clustStatType  {numStats}              = cell array of strings defining each cluster stat and correction method used (e.g.,  'size: vs biggest', 'size: stepDown', 'mass: multiVar FDR', 'mass: multivar Liptak')
%   sigClust_nodes {numStats}{numSigClust} = members cell array, each entry is list of ROI indices making up a significant cluster (based on sum cluster stat)
%   sigClust_pVal  {numStats}(numSigClust) = pvalue for each cluster in matched "node" list
%

%- return members of significant ROI clusters and the cluster p values
%sigUncor_nodes    = allmembers_true{1};  %- not cluster corrected
%iStat = 1;
%clustStatType{1}  = sprintf('%s: vs biggest',whichStat{iStat});
%sigClust_pVal{1}  = {};
%sigClust_nodes{1} = allmembers_true{1}(pCsum_true<=0.05){clustStatsInd_true{iStat}}


% sigroi_mCsum = allmembers_true{1}(pCsum_true<=0.05);  %- looking fo significant cluster stats at p<0.05 (DONT CHANGE this p threshold)
% sigroi_pCsum = pCsum_true(pCsum_true<=0.05);
% sigroi_mCsiz = allmembers_true{1}(pCsiz_true<=0.05);  %- looking fo significant cluster stats at p<0.05 (DONT CHANGE this p threshold)
% sigroi_pCsiz = pCsiz_true(pCsiz_true<=0.05);


% %- create summary string (good for command line output or plot labels)
% ROInotSigBoth = setxor([sigroi_mCsum{:}],[sigroi_mCsiz{:}]); if length(ROInotSigBoth)>0,strNotSame=',ROIsize~=sum'; else strNotSame=''; end
% %strOut=sprintf('\n Cluster (%d ROI, %d perm, sigThresh=%.2f, lumpPosNeg=%d): <%d sig size,%d sig sum>',numROI,numPerm,sigThresh,lumpPosNeg,length(sigroi_mCsiz),length(sigroi_mCsum));
% strOut=sprintf('\n Sig %s Clusters <Sum=%d,Size=%d%s>',strLump,length(sigroi_mCsum),length(sigroi_mCsiz),strNotSame);
% for iClus=1:length(sigroi_mCsum), strOut=sprintf('%s [%d] %d ROIs, p=%.04f;',strOut,iClus,length(sigroi_mCsum{iClus}),sigroi_pCsum(iClus)); end;
% if length(sigroi_mCsum)==0,  strOut=sprintf('%s NO SIG CLUSTERS',strOut); end;
% 
% if length(ROInotSigBoth)>0,
%     strOut2=sprintf('\n Size clusters:');
%     for iClus=1:length(sigroi_mCsiz), strOut2=sprintf('%s [%d] %d ROIs, p=%.04f;',strOut2,iClus,length(sigroi_mCsiz{iClus}),sigroi_pCsiz(iClus)); end;
%     fprintf('%s',strOut2);
% end
% 
% fprintf(' [%.1f s]',toc);
% if size(permdata,2)>1, disp(strOut); else disp(' '); end % no display if just a single perm... using this function to get uncorrected clusters in that case



%-- sanity check -- confirmed that emperic is symetric (prob not true if permuting over trials, but if permuting over subjects then it is)
% permdataOG = permdata;
% permdata                           = permdataOG;
% permdata(isnan(permdata))          = 0;
% permdata(     (permdata) < sigThresh) = 0; %-pos (no negs)
% [CsumMax_permPos,CsizMax_permPos,allmembers_permPos] = getSigROI_clusterCorrectTstatsJW2(permdata,adj);
%
% permdata                           = permdataOG;
% permdata(isnan(permdata))          = 0;
% permdata(    -(permdata) < sigThresh) = 0; %-neg (no pos)
% [CsumMax_permNeg,CsizMax_permNeg,allmembers_permNeg] = getSigROI_clusterCorrectTstatsJW2(permdata,adj);
%
% figure; subplot(311); hist(CsizMax_permPos,20); subplot(312); hist(CsizMax_permNeg,20);  subplot(313); hist(CsizMax_perm,20);
% figure; subplot(311); hist(CsumMax_permPos,20); subplot(312); hist(CsumMax_permNeg,20);  subplot(313); hist(CsumMax_perm,20);

