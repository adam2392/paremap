function [Hclust HsigClust NumPermOut StatStr] = perfClusterStats( subjXfreqXtime, pThreshCluster, UNCORR_CLUST );
%
%  perfClusterStatsSubj is used to do a cluster correction on across subject data, where each datapoint is the mean difference from one subject
%      i.e., randomly selecting which subjects will have all trial types ("correct" vs "incorect") swapped to see if across subject effect survives
%      (as opposed to the computationally challenging approach of permuting trial labels within subject to create a distribution  of mean difference for each subject
%
%
% INPUTS:
%  subjXfreqXtimeXfreq is mean DIFFERENCES or T-STAT: can be a time series or a spectrogram compairing two conditions
%         dim1: rows = subjects;    dim2: columns = freqs (=1 if just time series);    dim3: timepoints;
%
%  pThreshCluster: percentile used to binarize the data and define clusters;  actual cluster stat is always reported at 0.05
%
% OUTPUTS:
%  Hclust:     nan matrix with 1 for cluster locations (not multiple comparison corrected)
%  HsigClust:  same as Hclust, but only with clusters that survive multiple comparisons
%  NumPermOut: number of permutations used in this computation
%  TstatSpec:  true t-stat across the spectrogram


tstart=tic; tic;
drawnow; currentFig=gcf; %- save so the control can be passed back to the current fig after execution


%
%  I.   Compute the true mean difference
%  II.  Compute the distribution of mean differnces by permuting labels --- in this case that means randomly flipping the sign for each subject
%  III. Compute the per
%  IV.  Cluster correct the numPerm
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-  USER SELECTED OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SHOW_CLUSTER_FIG=0;  %- (0) defulat:  set to 1 to see intermediate steps.  Also can tinker with  separatePosNeg and collapsePosNeg to see additional plots of those variants
FIG_OFFSET = 200;


numPermMax = 1e6;
numPermMax = 2000;


%- OPTIONS: lump positive and negative data (just look for clusters of absolute differences), or find clusters separately for each.
separatePosNeg  = 1;     %- 0=dont separate (lump);  1=analyze each separately (efficient version)

%- MINIMUM CLUSTER SIZE: 
minClustSize    = 10;    %- (1,5,10) minimum cluster size for calculating cluster statistics.  for ROI's use 1, for Spectrograms use ~5 or 10; if pThresh=0.1 use 20; if pThresh=0.05 use 10?


%- PICK ONE OF THESE METHODS
GET_P_rankPerm  = 0;     %- straight-up ranking the permutation sums... most accurate p-value (valid for any distribution)
GET_P_zScrPerm  = 0;     %- z-score permuted values                                                     --- fastest method
GET_P_tTest     = 1;     %- t-test of actual values (and permuted t-statistics for cluster correction)  --- about 2x slower than rankPerm & zScore


if UNCORR_CLUST,
    numPermMax     = 1;  %- no reason to permute if just getting pValue from t-test
    GET_P_rankPerm = 0;  %- straight-up ranking the permutation sums... most accurate p-value (valid for any distribution)
    GET_P_zScrPerm = 0;  %- z-score permuted values                                                     --- fastest method
    GET_P_tTest    = 1;  %- t-test of actual values (and permuted t-statistics for cluster correction)  --- about 2x slower than rankPerm & zScore
    %minClustSize  = 1;  %- show all significant pixels, not just minClustSize 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Initailze parameters and check options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('pThreshCluster','var')
    pThreshCluster = 0.05;  %- standard is 0.05; reduce specificity by increasing
end

% OPTIONS:        'CLUSTER'     -- cluster correction
MULT_COMPARE_TYPE = 'CLUSTER';  %- FDR, SIMPLE_PERM, CLUSTER
multCompPthresh   = 0.05;       %- significant pValue for TRUE cluter relative to perm clusters.  ALWAYS keep at 0.05...


%- parameter check
if GET_P_rankPerm+GET_P_zScrPerm+GET_P_tTest~=1, fprintf('\nERROR: exxactly 1 p-value type should selected, defaulting to GET_P_rankPerm'); keyboard; GET_P_rankPerm=1; GET_P_zScrPerm=0; GET_P_tTest=0; end


%- output stats string that reports the method being used
if GET_P_rankPerm>0, strPval = 'pRank'; end
if GET_P_zScrPerm>0, strPval = 'zScor'; end
if GET_P_tTest>0,    strPval = 'tTest'; end

if separatePosNeg>0, strPN = 'SepPN'; else strPN = 'LumpPN'; end
%strPN = sprintf('%s%s',strPval,strPN);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- DEFINE THE PERMUTATIONS,  each subject gets a +1 or -1 (randomly assigned) for each permutation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get counts for each dimension
numSubj    = size(subjXfreqXtime,1);
numFreq    = size(subjXfreqXtime,2);
numTime    = size(subjXfreqXtime,3);
numPermAll = 2^numSubj;


%- EXACT (EXAUSTIVE) PERMUTATION TEST -- (1) should include all permutations including the "observed" effect and         (2) permutations are drawn without replacement (Smyth & Phipson 2010)
%-     RANDOM-SAMPLE PERMUTATION TEST -- (1) should include the observed effect to implicitly compute the EXACT p-value  (2) permutations are drawn without replacement (Smyth & Phipson 2010)

%- randomly pick a subset if "numPermMax" is less than "numPermAll"
if numPermAll>numPermMax,
    iPermUse   = sort([1 1+randperm(numPermAll-1,numPermMax)],'ascend'); %- keep "true" as first index
    numPermUse = numPermMax+1;
    NumPermOut = numPermMax;
else
    iPermUse   = [1:numPermAll];
    numPermUse = numPermAll;
    NumPermOut = numPermAll;
end

%- for each permutation (between 1 and numPermMax) assign each subject a +1 or -1.
%- follwoing method draws permutations WITHOUT REPLACEMENT by using a unique decimal-2-binary conversion for each possible permutation
permList = nan(numSubj,numPermUse); %- initialize matrix with nan; these should all be replaced by zeroes and ones and converted to -1 or 1
for iPerm=[1:numPermUse],
    permList(:,iPerm) = rem(floor((numPermAll-iPermUse(iPerm))*pow2(1-numSubj:0)),2);  %- algorithm pulled from the guts of "dec2bin"
end
permList=permList*2-1;


%- create output string that defines the type of stat computed here
StatStr = sprintf('CLUSTER: %s, pThr %.3f; %s; minClus %d; %d perm',strPval,pThreshCluster,strPN, minClustSize, NumPermOut); 
StatStr(find(StatStr=='_'))=' ';


%- create flattened version of spectrogram used in subsequent steps
subjXfreqXtimeFlat = reshape(subjXfreqXtime,[numSubj numFreq*numTime]);  % flatten time x freq for matrix manipulations 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- CALCULATE P-VALUES for EACH PERMUTATION: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GET_P_tTest==1,        %%-- Compute P-Values either by Ranking by fitting to a normal distribution
    
    effect_perm = nan(numPermUse,numFreq*numTime);
    onesMat     = ones(1,numFreq*numTime);
    
    if separatePosNeg>0, tailStr='right'; else tailStr='both'; end %- 'right' is 1-tailed test and will make very negative values have a p near 1.0, very positive near 0.0
    for iPerm=1:numPermUse, %- less perms than timexfreq points... so loop over perms
        [H,P,CI,STATS] = ttest(subjXfreqXtimeFlat.*(permList(:,iPerm)*onesMat),0,'dim',1,'tail',tailStr);  %- test whether mean is different from zero
        effect_perm(iPerm,:)    = STATS.tstat;
        effect_pValUse(iPerm,:) = P;
    end
    zTrans_perm      = effect_perm; %- used in ROI cluster correction machinery
    effect_perm      = reshape(effect_perm,[numPermUse numFreq numTime]);
    effect_pValUse   = reshape(effect_pValUse,[numPermUse numFreq numTime]);
    
    %- define a threshold value for the effect data --- different for z-score vs t-stat
    sigThresh = tinv(1-pThreshCluster/(separatePosNeg+1),numSubj-1); %- divide by 2 if separately analyzing positive and negative; degrees of freedom = numSubj-1

    
elseif GET_P_rankPerm,       %%-- Compute P-Values by Ranking empirical distribution (accurate for any distrubtion)
    
    effect_perm = (subjXfreqXtimeFlat'*permList)'/numSubj;   % result= PERMs x TIME_FREQ:   ~40x faster than loop of mean calculations
    
    %- compute a p-value by ranking the data and saying how many points are above or below:  [slowest step]
    effect_perm   = reshape(effect_perm,[numPermUse numFreq numTime]);
    effect_pValRank = nan(size(effect_perm));  % (numPerm x freq x time)
    for iFreq=1:numFreq,
        for iTime=1:numTime,
            if separatePosNeg,
                [uSort, ~,iRank] = unique(    (effect_perm(:,iFreq,iTime)),'sorted' );  %- raw value (not abs) to separate positive and negative effects
                effect_pValRank(:,iFreq,iTime) = (1-((iRank-0.5)/length(uSort)));     %- use length uSort (instead of numPerm) to deal with repeated values; p = 1-fraction of values that are smaller than this one
            else
                [uSort, ~,iRank] = unique( abs(effect_perm(:,iFreq,iTime)),'sorted' );  %- absolute value to convert to a 1-tailed test (input is difference, so deviation from zero is what we are looking for
                effect_pValRank(:,iFreq,iTime)   = (1-((iRank-1)/length(uSort)));       %- use length uSort (instead of numPerm) to deal with repeated values; p = 1-fraction of values that are smaller than this one
            end
        end
    end
    effect_pValUse   = effect_pValRank;
    
    
    
elseif GET_P_zScrPerm,    %%-- Compute P-Values either by z-scoring the empirical distribution (assumes normality)
    
    effect_perm = (subjXfreqXtimeFlat'*permList)'/numSubj;   % mean effect = PERMs x TIME_FREQ:   ~40x faster than loop of mean calculations
    
    %- compute a p-value assuming the distibution at each time point is normal: z-transform and use cdf of normal distribution   --- [3rd slowest step]
    zMu = mean(effect_perm(2:end,:),1);  %- dont use "true" effect in zMu and zSD calculation
    zSD = std(effect_perm(2:end,:),0,1);
    zTrans_perm = effect_perm -  (zMu' * ones(1,numPermUse))';
    zTrans_perm = zTrans_perm ./ (zSD' * ones(1,numPermUse))';
    
    if separatePosNeg, effect_pvalZscr =   (1-normcdf(   (zTrans_perm),0,1));      % get p-Values from raw z-score; used for separate positive and negative effects
    else               effect_pvalZscr = 2*(1-normcdf(abs(zTrans_perm),0,1)); end  % get p-values from abs(zscore); two-tailed test at 0.05 is actually 0.025 on each side, so multiply p of abs value by 2
    
    %- reshape to freq by time
    effect_pValUse = reshape(effect_pvalZscr,[numPermUse numFreq numTime]);
    
    
    %- define a threshold value for the effect data --- different for z-score vs t-stat
    sigThresh = norminv(1-pThreshCluster/(separatePosNeg+1));          %- divide by 2 if separately analyzing positive and negative
    
    
end

fprintf('<stats%.1f', toc(tstart)); tic;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Run the cluster correction algorithm and package the result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
USE_ZAR_METHOD = 0;
if USE_ZAR_METHOD,
    clear p_sig p_boot
    p_sig(1,:,:)    = effect_pValUse(1,:,:);      % psig  = (channels,frequency,time)
    p_boot(:,1,:,:) = effect_pValUse(2:end,:,:);  % pboot = (Nperm,chan,frequency,time)
    if separatePosNeg==0, clusSep=0; else clusSep=2; end
    
    [clusPvalPstat, clusPos, clusPos_sigstat, clusPos_sigsize] = getSigClustSpectrogram(p_sig, p_boot, pThreshCluster, multCompPthresh, clusSep, false, SHOW_CLUSTER_FIG*(FIG_OFFSET+clusSep));
    
    fprintf('|%.1fs>', toc); tic
    
    %- package outputs for return
    clusPos         = squeeze(clusPos);                 %- cluster positions
    clusPos_sigstat = squeeze(clusPos_sigstat);         %- cluster p-values
    clusPos_sigsize = squeeze(clusPos_sigsize);         %- cluster p-values
    Hclust          = clusPos./clusPos;                 %- true/false... thresholded cluster exists? (dimensions of p_sig)
    HsigClust       = clusPos_sigstat./clusPos_sigstat; %- true/false... significant cluster?        (dimensions of p_sig)
    %HsigClust       = (clusPos_sigstat+clusPos_sigsize)./(clusPos_sigstat+clusPos_sigsize); %- true/false... significant cluster?        (dimensions of p_sig)
else
    %- make fake unthresholded data for plots below
    clusPvalPstat   = squeeze(effect_pValUse(1,:,:));
    Hclust          = nan(numFreq,numTime);
    if separatePosNeg==0, Hclust(clusPvalPstat<=pThreshCluster) = 1;                                             %- div pThreshCluster by 2 if separatePosNeg = 0 (which should always be true)
    else                  Hclust(clusPvalPstat<=pThreshCluster/2 | clusPvalPstat>=(1-pThreshCluster/2)) = 1; end %- div pThreshCluster by 2 if separatePosNeg = 0 (which should always be true)
    HsigClust       = Hclust;
    clusPvalPstat   = 0.05*ones(sum(~isnan(Hclust(:))),1); %- just a list of number of significant nodes
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Run the cluster correction algorithm and package the result   --- ROI MACHINERY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[allClust_uncor, sigClust_stat, sigClust_ind, sigClust_pVal, strOut] = getSigClustA_mainROIs(truedata,permdata,adj,sigThresh,lumpPosNeg);
%   truedata (ROI x 1)
%   permdata (ROI x numPerm)
%   adj      (ROI x ROI)  - adjacency matric
%   sigThresh  (double)   - threshold value used to convert true and perm data into 0 and sig values
%   lumpPosNeg (0 or 1)   - 1=adjacent significant ROIs cluster irrespective of sign;  2=adjacent sig ROIs of same sign cluster

%%- USE ROI cluster machinery.
truedata = zTrans_perm(1,:)';     %- getSigClust expects nodes x perm
permdata = zTrans_perm(2:end,:)'; %- getSigClust expects nodes x perm
adj      = zeros(numFreq*numTime,numFreq*numTime);
for iTF=1:numFreq*numTime,
    
    thisFreq = mod(iTF-1,numFreq) + 1;
    thisTime = 1+floor((iTF-1)./numFreq);
    
    %- list of indicies that will be offset from iTF; respect spectrogram boundaries
    iC = 0; %- self connection
    
    %- 4-connected condition
    if thisFreq>1,       iC=[iC -1]; end       %- bottom
    if thisFreq<numFreq, iC=[iC +1]; end       %- top
    if thisTime>1,       iC=[iC -numFreq]; end %- left
    if thisTime<numTime, iC=[iC +numFreq]; end %- right
    
    %- also add these 4 corners for 8-connected condition
    USE_8CONN = 0;
    if USE_8CONN,
        if thisFreq>1&thisTime>1,             iC=[iC -1-numFreq]; end %- bottom-left
        if thisFreq>1&thisTime<numTime,       iC=[iC -1+numFreq]; end %- bottom-right
        if thisFreq<numFreq&thisTime>1,       iC=[iC +1-numFreq]; end %- top-left
        if thisFreq<numFreq&thisTime<numTime, iC=[iC +1+numFreq]; end %- top-right
    end
    
    iC = iC+iTF;
    adj(iTF,iC)=1; 
end
%figure(99); clf; imagesc(reshape(adj(105,:),[numFreq numTime])); set(gca,'ydir','normal'); %- sanity check that adjacency matrix is correct

%sigThresh   = 1.65;  %- (1.96) for 2-tailed 0.05;   (1.65) for 2-tailed 0.10
%sigThresh   = norminv(1-pThreshCluster/2);      %- use same threshold as old method... divide by 2 because lumpPosNeg should always be 0 (if 1, then dont div 2)
%sigThresh   = tinv(1-pThreshCluster/2,numSubj-1); %- use same threshold as old method... divide by 2 because lumpPosNeg should always be 0 (if 1, then dont div 2)
lumpPosNeg   = xor(separatePosNeg,1);  %- lumpPosNeg should = 0; always evaluate pos and negative separately 

[allClust_uncor, sigClust_stat, sigClust_ind, sigClust_pVal, strOut] = getSigClustA_main(truedata,permdata,adj,sigThresh,lumpPosNeg,minClustSize);


iStat = 1;
allNodes   = [allClust_uncor{:}];
sigNodes   = [allClust_uncor{sigClust_ind{iStat}}];
%sigNodes  = [allClust_uncor{4}]; %- sanity check
Hclust2    = nan(numFreq,numTime);
HsigClust2 = nan(numFreq,numTime);
Hclust2(allNodes)   =1;
HsigClust2(sigNodes)=1;


%- little figure for comparing old and new cluster methods
figure(1999);clf; set(gcf,'position',[795  700  1200  600]);
%- all clusters
subplot(231); imagesc(Hclust);     title(sprintf('OG all (%d clust; pThr %.3f)',length(find(clusPvalPstat)),pThreshCluster));    set(gca,'ydir','normal');
subplot(232); imagesc(Hclust2);    title(sprintf('new all (%d clust; %d minSize)',length(allClust_uncor),minClustSize));   set(gca,'ydir','normal');            
Hcomb=Hclust2; Hcomb(isnan(Hcomb))=0;Hcomb(find(Hclust==1))=Hcomb(find(Hclust==1))+1;
subplot(233); imagesc(Hcomb);      title('overlap');   set(gca,'ydir','normal');
%- sig clusters
subplot(234); imagesc(HsigClust);  title(sprintf('OG sig sum (%d clust)',length(find(clusPvalPstat<=0.05))));         set(gca,'ydir','normal');
subplot(235); imagesc(HsigClust2); title(sprintf('%s (%d clust)',sigClust_stat{iStat},length(sigClust_ind{iStat})));  set(gca,'ydir','normal');  
Hcomb=HsigClust2; Hcomb(isnan(Hcomb))=0;Hcomb(find(HsigClust==1))=Hcomb(find(HsigClust==1))+1;
subplot(236); imagesc(Hcomb);      title('overlap');    set(gca,'ydir','normal');


%adj(:,:)=0; for iTF=1:length(adj), adj(iTF,iTF)=1; end %- use diagnal to convert from clusters to single voxel tests

%- add a drawnow so any cluster stat figs are updated before we point back to the calling function's figure
drawnow; 
pause(1);
%jwSaveFigs([1999 80]); %- fig 80 


HsigClust = HsigClust2;
figure(currentFig);
pause(1);


