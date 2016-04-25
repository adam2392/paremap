
%
%  perfClusterStatsTest_AL is used to do a cluster correction on within
%      subject data, where each datapoint is the mean of one specific trial
%      (e.g. BRICK, CLOCK,etc. probewords) swapped to see if each word
%      pairing's power effect survives 
%      (as opposed to the computationally challenging approach of 
%       permuting trial labels within subject to create a distribution  
%       of mean difference for each subject
%
clc;
clear all;
subj = 'NIH034';

addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation');

% start time
tstart=tic; tic;
drawnow; currentFig=gcf; %- save so the control can be passed back to the current fig after execution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/');

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

%%- GET CORRECT EVENTS ONLY
% POST MODIFY EVENTS based on fields we want (e.g. is it correct or not)?
correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);

%%- Check all target words == vocalized words
% vocalizedWords = {events.vocalizedWord};
% targetWords = {events.targetWord};

%?? questions
% 1. What is Uncorr_clust

%%- PROCESS
%
%  I.   Compute the true mean difference
%  II.  Compute the distribution of mean differnces by permuting labels 
%   --- in this case that means randomly flipping the sign for each EVENT
%   TRIGGER
%  III. Compute the per
%  IV.  Cluster correct the numPerm
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-  USER SELECTED OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SHOW_CLUSTER_FIG=0;  %- (0) defulat:  set to 1 to see intermediate steps.  Also can tinker with  separatePosNeg and collapsePosNeg to see additional plots of those variants
FIG_OFFSET = 200;

% start with 500, then 2000
numPermMax = 500; % maximum # of permutations

%- MINIMUM CLUSTER SIZE: 
minClustSize    = 5;    %- (1,5,10) minimum cluster size for calculating cluster statistics.  for ROI's use 1, for Spectrograms use ~5 or 10; if pThresh=0.1 use 20; if pThresh=0.05 use 10?

%- PICK ONE OF THESE METHODS
GET_P_rankPerm  = 0;     %- straight-up ranking the permutation sums... most accurate p-value (valid for any distribution)
GET_P_zScrPerm  = 0;     %- z-score permuted values                                                     --- fastest method
GET_P_tTest     = 1;     %- t-test of actual values (and permuted t-statistics for cluster correction)  --- about 2x slower than rankPerm & zScore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Initailze parameters and check options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('pThreshCluster','var')
    pThreshCluster = 0.05;  %- standard is 0.05; reduce specificity by increasing
end
%- output stats string that reports the method being used
if GET_P_rankPerm>0, strPval = 'pRank'; end
if GET_P_zScrPerm>0, strPval = 'zScor'; end
if GET_P_tTest>0,    strPval = 'tTest'; end

% OPTIONS:        'CLUSTER'     -- cluster correction
MULT_COMPARE_TYPE = 'CLUSTER';  %- FDR, SIMPLE_PERM, CLUSTER
multCompPthresh   = 0.05;       %- significant pValue for TRUE cluter relative to perm clusters.  ALWAYS keep at 0.05...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- DEFINE THE PERMUTATIONS,  each trial within subject gets a +1 or -1 (randomly assigned) for each permutation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/freq_probeToVocal_100msbinned/1_G1-global_probeToVocal.mat';
data = load(file);
data = data.data;
fields = fieldnames(data);
expXfreqXtime = data.powerMatZ;

%%- get counts for each dimension
numExps    = size(expXfreqXtime, 1);    % # of events (within subjects)
numChans   = 96; % there are 96 channels, but we will try to save them separately
numFreq    = size(expXfreqXtime,2);     % # of frequency indices
numTime    = size(expXfreqXtime,3);     % # of time indices
numPermAll = 2^numChans; %factorial(numChans);        % 2^numSubjs if doing across subjects

%- EXACT (EXAUSTIVE) PERMUTATION TEST -- 
%   (1) should include all permutations including the "observed" effect and         
%   (2) permutations are drawn without replacement (Smyth & Phipson 2010)
%- RANDOM-SAMPLE PERMUTATION TEST -- 
%   (1) should include the observed effect to implicitly compute the EXACT p-value  
%   (2) permutations are drawn without replacement (Smyth & Phipson 2010)
%%%%% ?? NumPermOut vs. numPermUse ?? %%%%%%
%- randomly pick a subset if "numPermMax" is less than "numPermAll"
if numPermAll>numPermMax, % RANDOM-SAMPLE PERMUTATION TEST
    if numPermAll > 2^53
        numPermAll = 2^52;
    end
    iPermUse   = sort([1 1+randperm(numPermAll-1,numPermMax)],'ascend'); %- keep "true" as first index
    numPermUse = numPermMax+1;
    NumPermOut = numPermMax;
else                      % EXHAUSTIVE PERMUTATION TEST
    iPermUse   = [1:numPermAll];    % the index of permutation
    numPermUse = numPermAll;        % # of permutations to use
    NumPermOut = numPermAll;        % # of permutations used
end

%%%%%%%%%%%%%%%%% What should I implement here? Grab two unique event
%%%%%%%%%%%%%%%%% triggers and relable?
%- for each permutation (between 1 and numPermMax) assign each trial a +1 or -1.
%- following method draws permutations WITHOUT REPLACEMENT 
%   by using a unique decimal-2-binary conversion for each possible permutation
permList = nan(numExps,numPermUse); %- initialize matrix with nan; these should all be replaced by zeroes and ones and converted to -1 or 1
for iPerm=[1:numPermUse],
    permList(:,iPerm) = rem(floor((numPermAll-iPermUse(iPerm))*pow2(1-numChans:0)),2);  %- algorithm pulled from the guts of "dec2bin"
end
permList=permList*2-1;

%- OPTIONS: lump positive and negative data (just look for clusters of absolute differences), or find clusters separately for each.
separatePosNeg  = 1;     %- 0=dont separate (lump);  1=analyze each separately (efficient version)

if separatePosNeg>0, strPN = 'SepPN'; else strPN = 'LumpPN'; end
%- create output string that defines the type of stat computed here
StatStr = sprintf('CLUSTER: %s, pThr %.3f; %s; minClus %d; %d perm',strPval,pThreshCluster,strPN, minClustSize, NumPermOut); 
StatStr(find(StatStr=='_'))=' ';