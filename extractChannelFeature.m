%%%%% Function: extractChannelFeature.m
%%%%% Description: Compute a distance metric from the condensed_data
%%%%% directory saved output. The data is in time binned and freq. binned
%%%%% matrix.
%%%%% 
%%%%% Input: 
%%%%% channel_num - the channel number we want to extract (e.g.
%%%%% 1-96)
%%%%% probeword - 
%%%%% USEPROBETOVOCAL - either 1, or 0. If we want to extract the
%%%%% probetolocal data struct instead of the whole data struct
%%%%% 
%%%%% Output:
%%%%% datapoints - vector of datapoints (feature) for this channel
%%%%%
function [features] = extractChannelFeature(channel_num, USEPROBETOVOCAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data that can plot evoked and spectrogram
% data.uniqueTrigType = uniqueTrigType; % store all the unique trigger types
% data.trigType = trigType;             % store the trigger type per event
% data.wavesSft = wavesSft;             % store teh eeg wave forms
% data.waveT = waveT;                   % time points for the wave form
% data.powerMatZ = powerMatZ;           % z-scored power
% data.waveletFreqs = waveletFreqs;     % wavelet frequencies used
% data.chanNum = channel_num;           % store the corresponding channel number
% data.chanStr = chanStr;               % the string name of the channel
% data.freqBandYtick = freqBandYticks;
% data.freqBandYlabel = freqBandYtickLabels;

%% Load the Data From Dir
if ~exist('data')
    channel_num = '1';
    chanStr = 'G1-global';
    dataDir = 'condensed_data/';
    filename = strcat(dataDir, channel_num, '_', chanStr);
end

data = load(filename);
data = data.data;

%% SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];

eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

% only get the correct events
correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);
clear correctIndices 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load Data and Set Triggers     ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Load in data parameters we need to do statistical testing
powerMatZ           = data.powerMatZ;
waveT               = data.waveT;
trigType            = data.trigType;
uniqueTrigType      = data.uniqueTrigType;
chanNum             = data.chanNum;
chanStr             = data.chanStr;
freqBandYticks      = data.freqBandYtick;
freqBandYtickLabels = data.freqBandYlabel;
waveletFreqs        = data.waveletFreqs;

%% GET TRIGGERS WE WANT
%%- As of 03/04/16: get each trigger word by itself, then word A with B and
%%C and then word B with C
% loop through each trigger and create events struct to pass to plotting
TRIGGER_TYPES   = {'BRICK','CLOCK','GLASS','JUICE','PANTS'}; %-fixation and blockStart not ready for physio analysis
sampEventsMeta = events;  % includes assocaited + and *
probeWords = {sampEventsMeta.probeWord};
targetWords = {sampEventsMeta.targetWord};

%%- Loop through each probeword
for i=1:length(TRIGGER_TYPES)
    THIS_TRIGGER = TRIGGER_TYPES{i}; % set the current probeword
    
%     THIS_TRIGGER = probeword;
    %%- 01: GET TRIGGER INDICES WE WANT
    switch THIS_TRIGGER,
        %%- For each probeword:
        % - find events with that probeword
        % - get the unique targetwords for that event
        case 'BRICK'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'BRICK PROBE';

            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
        case 'CLOCK'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'CLOCK PROBE';
            
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
       case 'JUICE'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'JUICE PROBE';
            
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
       case 'PANTS'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'PANTS PROBE';
            
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
       case 'GLASS'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'GLASS PROBE';
            
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
        otherwise
            error('no event trigger selected');
    end
    
    %%- 02: GO THROUGH EACH TARGETWORD FOR THIS PROBEWORD (THISTRIGGER)
    for j=1:length(targets) % loop through each unique trigger for a specific probeword
        % find event indices for this trigger matched with a specific
        % targetword
        targetWord = targets{j};
        eventInd = find(strcmp({sampEventsMeta.probeWord},THIS_TRIGGER) & strcmp({sampEventsMeta.targetWord},targetWord));
    
        %%%%%% For each matrix, get a data point measure, store it in a
        %%%%%% vector. Then for each probe_word 
        %%- pass in events, powerMatZ, indices of events we want
        %%- Spectrograms of each one
        thisPowMat = powerMatZ(eventInd,:,:);
        powPlot = mean(thisPowMat, 1); % average across events
        powPlot = squeeze(powPlot);
        
        % add rows/fields to the features struct
        
        
    end
end

%%- Save data struct, and/or return it
% now that we have a struct of features for each probeword/target_word
% (match)
% compute centroid of each match and compute distances from each centroid

