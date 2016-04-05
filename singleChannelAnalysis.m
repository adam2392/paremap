
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
clc;
%% SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];

%% EXTRACTION OPTIONS
TRIGGER_TYPES   = {'BRICK','CLOCK','GLASS','JUICE','PANTS'}; %-fixation and blockStart not ready for physio analysis
THIS_TRIGGER    = TRIGGER_TYPES{1};   %%%%%%%%%% SELECT TRIGGER TYPE FOR EXTRACTION AND ANALYSIS
%% PLOTTING OPTIONS
HIDE_FIGURES    = 0;
USE_CHAN_SUBSET = 0; %0=all channels (not the subset); >=1 means process than many of the subset
FIG_OFFSET = 0;

%%- PLOT PARAMETERS
figFontAx       = 18;    if ispc, figFontAx = figFontAx-5; end
if ~exist('FIG_OFFSET','var'), FIG_OFFSET = 0; end %- default to 0, but if called as function then allow that value to persist

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

%%-array of frequency bands
freqBandAr(1).name    = 'delta';
freqBandAr(1).rangeF  = [2 4];          %[2 4]
freqBandAr(2).name    = 'theta';
freqBandAr(2).rangeF  = [4 8];          %[4 8]
freqBandAr(3).name    = 'alpha';
freqBandAr(3).rangeF  = [8 16];         %[8 12]
freqBandAr(4).name    = 'beta';
freqBandAr(4).rangeF  = [16 32];        %[12 30]
freqBandAr(5).name    = 'low gamma';
freqBandAr(5).rangeF  = [32 80];        %[30 70]
freqBandAr(6).name    = 'high gamma';
freqBandAr(6).rangeF  = [80 160];       %[70 150]
freqBandAr(7).name    = 'HFO';
freqBandAr(7).rangeF  = [160 400];      %[150 400]

% set the frequency bands to certain ranges for plotting
for iFB=1:length(freqBandAr),
    freqBandAr(iFB).centerF = mean(freqBandAr(iFB).rangeF);
    %freqBandAr(iFB).label   = sprintf('%s-%.0fHz', freqBandAr(iFB).name(1:[ min( [length(freqBandAr(iFB).name), 6] )]), freqBandAr(iFB).centerF);
    freqBandAr(iFB).label   = sprintf('%s [%.0f-%.0f Hz]', freqBandAr(iFB).name, freqBandAr(iFB).rangeF);
end

freqBandYticks  = unique([freqBandAr(1:7).rangeF]);
for iFB=1:length(freqBandYticks), freqBandYtickLabels{iFB} = sprintf('%.0f Hz', freqBandYticks(iFB)); end

% clear variables to develop easier...
clear chanListUse chanStrUse
clear chanFile chanNums chanTags iChanList jackSheet iFB 
clear docsDir eegRootDir eegRootDirHome eegRootDirWork talDir behDir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load Data From Preprocessed Dir       ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condensedDataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/time_freq_binned/';
ext = '*.mat';
files = dir(fullfile(condensedDataDir, ext));
files = {files.name};

% for fileI=1:length(files),
    

% 82 and 33 max min
file = fullfile(condensedDataDir,'82_LP2-global.mat');
data = load(file);
data = data.data;

%% GET TRIGGERS WE WANT
% loop through each trigger and create events struct to pass to plotting
TRIGGER_TYPES   = {'BRICK','CLOCK', 'JUICE', 'PANTS', 'GLASS'}; %-fixation and blockStart not ready for physio analysis
TRIGGER_TYPES   = {'BRICK','GLASS'};
sampEventsMeta = events;  % includes assocaited + and *
probeWords = {sampEventsMeta.probeWord};
targetWords = {sampEventsMeta.targetWord};

close all
figMeta = figure(FIG_OFFSET+1); set(gcf, 'color', 'w')  % set color background to white
figSpect = figure(FIG_OFFSET+2); set(gcf, 'color', 'w')  % set color background to white
figAnova = figure(FIG_OFFSET+3); set(gcf, 'color', 'w')  % set color background to white

anovaPowMat = {}; % create a cell array to store all powerMatrices for certain trigger
nout = {};
xout = {};
%%- Loop through each probeword
for i=1:length(TRIGGER_TYPES)
    THIS_TRIGGER = TRIGGER_TYPES{i}; % set the current probeword
    
    %%- 01: GET TRIGGER INDICES WE WANT
    switch THIS_TRIGGER,
        %%- For each probeword:
        % - find events with that probeword
        % - get the unique targetwords for that event
        case 'BRICK'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'BRICK PROBE';
            
            tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
        case 'CLOCK'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'CLOCK PROBE';
            
            tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
       case 'JUICE'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'JUICE PROBE';
            
            tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
       case 'PANTS'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'PANTS PROBE';
            
            tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
       case 'GLASS'
            tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'GLASS PROBE';
            
            tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
            %%- get all the unique targetwords for BRICK probeword
            targets = unique({tempevents.targetWord});
        otherwise
            error('no event trigger selected');
    end
    
%     tempevents
%     length(tempInd)
    metaEvents = events(tempInd);
    
    %%- 02: Plot Meta Data For Each Trigger
    subplotIndex = i;
    NUMSUBPLOTS = length(TRIGGER_TYPES);
    n = plotEventMetaData(metaEvents, subplotIndex, figMeta.Number, NUMSUBPLOTS, THIS_TRIGGER);
    nout{i} = n;
%     xout{i} = x;
    
    %%- 03: Plot Spectrogram For Each Trigger
    powerMatZ = data.powerMatZ(tempInd,:,:);
    size(powerMatZ)
    triggers = data.trigType;
    % create time vector that is binned and still centered at 0
    binnedWaveT = [1:size(powerMatZ,3)] - round(data.timeZero);
    titleStr = sprintf('mean power: chan %s, %d events', data.chanStr, size(powerMatZ,1));
    
    plotEventSpectrogram(powerMatZ, binnedWaveT, freqBandAr, ...
        titleStr, subplotIndex, figSpect.Number, NUMSUBPLOTS, THIS_TRIGGER)

    %%- Create groups for ANOVA
    anovaPowMat{i} = data.powerMatZ(tempInd,:,:);    
end

%%- 04: Plot ANOVA For each Trigger
titleStr = strcat(sprintf('ANOVA: chan %s, COMPARING: ', data.chanStr), TRIGGER_TYPES{:});    
plotEventANOVA(anovaPowMat, binnedWaveT, freqBandAr, titleStr, ...
    figAnova.Number)

% c2 = combnk(TRIGGER_TYPES,2)
% test = zeros(size(c2))
for i=1:length(TRIGGER_TYPES)
    for j=i:length(TRIGGER_TYPES)
%         sprintf('%d and %d', i, j)
        if i~=j
            
            sprintf('Comparing: %s and %s', TRIGGER_TYPES{i}, TRIGGER_TYPES{j})
            [h, p] = kstest2(nout{i}, nout{j})
        end
    end
end
