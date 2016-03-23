%%%%% only for use if saved anova matrix in a separate .mat file

clear
clc
% close all

%% SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];
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

eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

% Either go through all the sessions, or a specific session
if sessNum == -1 | length(sessNum)>1, % all sessions
    disp('STEP 1: Going through all sessions')
    session = 'Meta Session [all]';
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/');
    sessStr = '[all]';
else                                  % one session
    disp('STEP 1: Going through one session')
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/', session);
    sessStr = sprintf('[%d]',sessNum);
end

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

% only get the correct events
correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);
clear correctIndices 

% anovaDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/anova/';
anovaDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/brickclock_glasspants/';

ext = '*.mat';
files = dir(strcat(anovaDir, ext));
files = {files.name};

%%% 
file = strcat(anovaDir, files{1});
data = load(file);
data = data.data;


%%- Initialize matrices/vectors for speed
anovaMats = zeros(length(files), size(data.anovaMat,1), size(data.anovaMat,2));
timeFreqImpactSpect = zeros(7,size(data.anovaMat,2));
freqImpact = zeros(7,1);
timeImpact = zeros(1,size(data.anovaMat,2));
size(anovaMats)

timeBin = 35;
timeImpact = zeros(1,timeBin);
timeFreqImpactSpect = zeros(7,timeBin);

responseTimes = data.responseTimes;

% loop through each channel file
for i=1:length(files)
    file = strcat(anovaDir, files{i});
    data = load(file);
%     data = data.anovaData;
    data = data.data;

    %%- Create spectrogram of p-values and convert to sig/non-significant
    spect = data.anovaMat;
    
    %%- cut down anova Mat
    spect = spect(:, 1:timeBin);
    
    spect(spect>0.05) = 1;
    spect(spect<=0.05) = 0.5;
    spect(spect==1) = 0;
    spect(spect==0.5) = 1;
    
    %%- Analyze impact of each time/freq bin
    freqImpact = freqImpact + sum(spect,2);
    timeImpact = timeImpact + sum(spect,1);
    
    timeFreqImpactSpect = timeFreqImpactSpect + spect;
    % concatenate into anova matrix
    anovaMats(i,:,:) = data.anovaMat;
    
    % plot spectrogram
%     spect = data.anovaMat;
%     spect(spect>0.05) = 1;
%     figure;
%     imagesc(spect);
%     colorbar;
end
freqImpact
timeImpact
figure
bar(freqImpact)
title('Impact of Different Frequency Bands')
xlabel('frequency bands from delta -> HFO')
% xlabel([freqBandAr.label])
figure
subplot(211)
bar(timeImpact)
title('impact of time bins')
xlabel('time bins')
subplot(212)
hist(responseTimes)
title('response times')
xlabel('response times')

%%- time freq impact spectrogram
figure;
imagesc(timeFreqImpactSpect/96)
xlabel('time bins')
ylabel('frequency bins')
title('Proportion that was significant over the 96 channels')
colormap(jet)
colorbar

% figure;
% spect = squeeze(mean(anovaMats, 1));
% imagesc(spect)

%%- Find most significant channel
minp = 5;
maxp = 0;
channel_num = 0;
maxchannel_num =0;
minchannel_num = 0;

for i=1:length(files)
    data = anovaMats(i,:,:);
    
    spect(spect>0.05) = 1;
    data(data<=0.05) = 0.5;
    data(data==1) = 0;
    data(data==0.5) = 1;
    
    summedData = sum(data(:))/96;
    if(summedData > maxp)
        maxp = summedData;
        maxchannel_num = i;
    end
    if(summedData < minp)
        minp = summedData;
        minchannel_num = i;
    end
%     if(min(data(:)) < minp)
%         minp = min(data(:));
%         channel_num = i;
%     end
end
% minp
% channel_num

minp
minchannel_num
maxp
maxchannel_num

%% PLOT SOME DATA ABOUT THE MOST IMPACTFUL CHANNEL 
%- Spectrogram of the two different triggers
%- 
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
%%------------------ STEP 2: Load Data From Preprocessed Dir       ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condensedDataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/freq_probeToVocal_100msbinned/';
ext = '*.mat';
files = dir(fullfile(condensedDataDir, ext));
files = {files.name};
file = strcat(condensedDataDir, files{1});
data = load(file);
data = data.data;

channelWeWant = '88';
filename = [strfind(files, channelWeWant)]
% filename = cellfun(strfind(files, channelWeWant))
% filename = cellfun(@(s) strfind(channelWeWant, s), files)

% 82 and 33 max min
file = fullfile(condensedDataDir, strcat(channelWeWant, '*.mat'));
data = load(file);
data = data.data;


% loop through each channel file
for filei=1:length(files)
    
    file = strcat(anovaDir, files{filei});
    data = load(file);
    data = data.data;
    
    disp(['on file: ', file])
    %% EXTRACT DAT THAT WE WANT
    thisChan = data.chanNum;
    thisChanStr = data.chanStr;
    powerMatZ = data.powerMatZ;
    trigType = data.trigType;
    responseTimes = data.responseTime;
    
    %% RUN ANOVA 
    %%- Data is already time and frequency binned
    % squeeze spectrogram into eventsXfreqsXtime if necessary
    spectMat = squeeze(powerMatZ);
    size(spectMat);
    
    %% GET TRIGGERS WE WANT
    sampEventsMeta = events;  % includes assocaited + and *
    probeWords = {sampEventsMeta.probeWord};
    targetWords = {sampEventsMeta.targetWord};
    TRIGGER_TYPES   = {'BRICK','CLOCK','GLASS','JUICE','PANTS'}; %-fixation and blockStart not ready for physio analysis

    anovaPowMat = {}; % create a cell array to store all powerMatrices for certain trigger
    anovaGroups = {}; % create cell array to hold the groups string
    index = 0;
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
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'BRICK PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
            case 'CLOCK'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'CLOCK PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'JUICE'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'JUICE PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'PANTS'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'PANTS PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'GLASS'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'GLASS PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
            otherwise
                error('no event trigger selected');
        end %end of switch
        
        %%- 02: GO THROUGH EACH TARGETWORD FOR THIS PROBEWORD (THISTRIGGER)
        for j=1:length(targets) % loop through each unique trigger for a specific probeword
            % find event indices for this trigger matched with a specific
            % targetword
            targetWord = targets{j};
            eventInd = find(strcmp({sampEventsMeta.probeWord},THIS_TRIGGER) & strcmp({sampEventsMeta.targetWord},targetWord));
            
            index = index + 1;
            %%- Create groups for ANOVA
            anovaGroups{index} = strcat(THIS_TRIGGER, '_', targetWord);
            anovaPowMat{index} = spectMat(eventInd,:,:);
        end
        
        %%- Create groups for ANOVA
%         anovaPowMat{i} = spectMat(tempInd,:,:);
    end
end