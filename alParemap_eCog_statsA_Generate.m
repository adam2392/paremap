%function alParemap_eCog_statsA_Generate( subj, sessNum, THIS_TRIGGER, THIS_REF_TYPE, FIG_OFFSET )
%
%         -- select subset of behavioral events
%         -- get waveform, power, phase for a time window around the event of interest
%

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

REF_TYPES       = {'noreref', 'bipolar', 'global'};
THIS_REF_TYPE   = REF_TYPES{3}; % (1) noreref, (2) bipolar, (3) laplacian

%% PLOTTING OPTIONS
HIDE_FIGURES    = 0;
USE_CHAN_SUBSET = 1; %0=all channels (not the subset); >=1 means process than many of the subset
FIG_OFFSET = 0;

%%- PLOT PARAMETERS
figFontAx       = 18;    if ispc, figFontAx = figFontAx-5; end
if ~exist('FIG_OFFSET','var'), FIG_OFFSET = 0; end %- default to 0, but if called as function then allow that value to persist

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

%% FILTERING OPTIONS
BP_FILTER_RAW                 = 1;  %-0 or 1: apply a bandpass filter to the raw traces (1-499 hz)
PROCESS_CHANNELS_SEQUENTIALLY = 1;  %0 or 1:  0 means extract all at once, 1 means sequentially

%% Testing
% TRIGGER_TYPES = TRIGGER_TYPES{1}; % testing

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

subjDir = fullfileEEG(eegRootDir,subj); % directory to subject (e.g. NIH034)
docsDir = fullfileEEG(subjDir,'docs');  % directory to the docs (electordes.m, tagNames.txt, etc.)
talDir  = fullfileEEG(subjDir,'tal');
defaultEEGfile = fullfileEEG('/Volumes/Shares/FRNU/data/eeg/',subj,'/eeg.reref/');  % default event eegfile fields point here... switch to local before loading

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------STEP 2: Create Channel List          -------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%- STEP 2: Manipulate variables for printing and display %%%%%%%%%%%%%%%%%
%%- Get 1. # of channels to use, 
%%-     2. list of channels
%%-     3. list of channel names
%%- select all channels, or part of the subset of channels

jackSheet = fullfileEEG(docsDir, 'jacksheetMaster.txt');
[chanNums chanTags] = textread(jackSheet,'%d%s%*s');

%%% always look at all electrodes... worry about "good" and "bad" later (bad means inter-ictal activity or seizure activity)
%- three referencing options:  noreref (should manually subtract reference channel), reref bioploar, and reref laplacian
chanStr = {};   % cell for all the channel names
chanFile = 0;   % file for the channels (e.g. ~/NIH034/tal/leads.txt) 
chanList = [];  % list of the channels (e.g. 1-96)
iChanListSub = []; % list of the subset of channels we want to analyze (e.g. [48 1])

switch THIS_REF_TYPE
    case 'noreref'  
    case 'bipolar'
    case 'global' % look at global electrodes
        fprintf('STEP 1: Using Global referencing\n');
        chanFile      = [talDir '/leads.txt'];
        chanList      = textread(chanFile,'%d'); % read in the list of channels nums
        
        % set the names for each channel
        for iChan=1:size(chanList,1),
            chanStr{iChan} = sprintf('%s-global', chanTags{find(chanNums==chanList(iChan))} );
        end
        eventEEGpath  = '/eeg.reref/';
        
        iChanListSub  = [46 1];            %G1, G2, LF1, AST1,
    otherwise
        fprintf('Error, no referencing scheme selected');
end
%%- subset channel I want to extract
iChanListSub  = [46, 1];
%%- select all channels, or part of the subset of channels
if USE_CHAN_SUBSET==0,
    iChanList = 1:size(chanList,1);  %all possible channels
else
%     iChanList = iChanListSub(1:min([length(iChanListSub) USE_CHAN_SUBSET]));    %select subset of channels
    iChanList = iChanListSub;
end

% what is this doing here?
chanListUse = [];  chanStrUse = {};
for iChan=iChanList,
    chanListUse(end+1,:) = chanList(iChan,:);
    chanStrUse{end+1}    = chanStr{iChan};
end

% reset variables and create list of channels and their corresponding names
chanList = chanListUse;
chanStr  = chanStrUse;
numChannels = size(chanList,1);

% get local copy of eegfile
% localeegfile = regexprep(eventTrigger(iEvent).eegfile,defaultEEGfile,fullfileEEG(subjDir,eventEEGpath))

% print statements for debugging and process checking
fprintf('\n');
fprintf('STEP 2 -- %d channels to process for %s : %s', numChannels, subj, session);
fprintf('\n');
disp('Variables to use here are:')
disp('chanList, chanStr, numChannels')
fprintf('\n');

% clear variables to develop easier...
clear chanListUse chanStrUse
clear chanFile chanNums chanTags iChanList jackSheet iFB 
clear docsDir eegRootDir eegRootDirHome eegRootDirWork talDir behDir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 3: CREATE EVENT TRIGGERS                  ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trigType = [];
eventsTrig = [];

% loop through each trigger and create events struct to pass to plotting
for i=1:length(TRIGGER_TYPES)
    THIS_TRIGGER = TRIGGER_TYPES{i};
%     THIS_TRIGGER = TRIGGER_TYPES;
    
    % get the events for a specific TRIGGER WORD
    %%- Only looking at probe words right now
    sampEventsMeta = events;  % includes assocaited + and *
    probeWords = {sampEventsMeta.probeWord};
    targetWords = {sampEventsMeta.targetWord};
    
    % events for specific TRIGGER
    sampEventsTrig = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
    
    % # of trigger types (e.g. words)
%     metaZeroMS = [eventsMeta.mstime];       % get metaZeroMS of each event
    eventsTrig = [eventsTrig, sampEventsTrig];            % get the events for that specific trigger
    trigWord = {sampEventsTrig.probeWord};  % get the trigger word for this sample
    trigType = [trigType, trigWord];      % get the trigger types
    
    switch THIS_TRIGGER,
        case 'BRICK'
            brickevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'BRICK PROBE';
            eventsTriggerXlim = [-2.25 5.25];
%             eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
        case 'CLOCK'
            clockevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'CLOCK PROBE';
            eventsTriggerXlim = [-2.25 5.25];
%             eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
        case 'GLASS'
            glassevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'GLASS PROBE';
            eventsTriggerXlim = [-2.25 5.25];
%             eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
        case 'JUICE'
            juiceevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'JUICE PROBE';
            eventsTriggerXlim = [-2.25 5.25];
%             eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
        case 'PANTS'
            pantsevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
            disp(['Looking at trigger: ', THIS_TRIGGER]);
            metaYstr = 'PANTS PROBE';
            eventsTriggerXlim = [-2.25 5.25];
%             eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
        otherwise
            error('no event trigger selected');
    end 
    
    disp(['Looking at ', metaYstr, ' in line ~196!'])
end

clear metaYstr trigWord probeWords sampEventsMeta sampEventsTrig
clear THIS_REF_TYPE THIS_TRIGGER
% print statements for debugging and process checking
fprintf('\n');
fprintf('STEP 3 -- %d events to process for %s : %s', length(eventsTrig), subj, session);
fprintf('\n');
disp('Variables to use here are:')
disp('eventsTrig, trigType, eventsTriggerXlim')
disp(['Variables to use here if looking at probe words are:'])
disp(['brickevents, glassevents, pantsevents, juicevents, clockevents'])

%%- GET CORRECT EVENTS ONLY
% POST MODIFY EVENTS based on fields we want (e.g. is it correct or not)?
correctIndices = find([events.isCorrect]==1);
trigType = trigType(correctIndices);
eventsTrig = eventsTrig(correctIndices);
events = events(correctIndices);

%%- Now get each unique word pairings A/B, A/C, A/D, B/C, B/D
% loop through each trigger type
currentUniqueTrigType = unique(trigType);
for itrig = 1:length(currentUniqueTrigType)
    trigger = currentUniqueTrigType(itrig); % current trigger
     
    matchTriggers = find(strcmp({eventsTrig.targetWord},trigger));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 4: DATA INPUT TO GETE_MS, MULTIPHASEVEC3 AND ZSCORE
%%------------------      AND SET UP POWER, POWERZ AND PHASE MATRICS        ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Input to gete_ms
%%- Dependent only on eventsTriggerXlim: These stay the same regardless of how we process events
eventTrigger = eventsTrig;
eventOffsetMS   = eventsTriggerXlim(1)*1000;      % positive = after event time; negative = before event time
eventDurationMS = diff(eventsTriggerXlim)*1000;   % duration includes offset (i.e., if offset -500 and duration 1000, only 500 ms post event will be prsented)

OffsetMS        = eventOffsetMS;     % positive = after event time; negative = before event time
DurationMS      = eventDurationMS;   % duration includes offset (i.e., if offset -500 and duration 1000, only 500 ms post event will be prsented)
BufferMS        = 1000;              % grab excess data before/after event window so filters don't have edge effect
resampledrate   = 1000;              % don't resample... keep the 1kHz native sample rate

%%- apply a bandpass filter raw data? (i.e. pre-filter the wave?)
if BP_FILTER_RAW==1,
    preFiltFreq      = [1 499];   %[1 499] [2 250]
    preFiltType      = 'bandpass';
    preFiltOrder     = 2;
    preFiltStr       = sprintf('%s filter raw; %.1f - %.1f Hz',preFiltType,preFiltFreq);
    preFiltStrShort  = '_BPfilt';
    FIG_OFFSET       = FIG_OFFSET+100+round(preFiltFreq(2));  %keep this empty to avoid any filtering of the raw data
else
    preFiltFreq      = []; %keep this empty to avoid any filtering of the raw data
    preFiltType      = 'stop';
    preFiltOrder     = 1;
    preFiltStr       = 'Unfiltered raw traces';
    preFiltStrShort  = '_noFilt';
end

%%- NEEDED FOR EVERY EVENTS
%%- remap event pointer from default (server) to local copy of the EEG data
for iEvent=1:length(eventTrigger),
    eventTrigger(iEvent).eegfile = regexprep(eventTrigger(iEvent).eegfile,defaultEEGfile,fullfileEEG(subjDir,eventEEGpath));

end

%%- gets the range of frequencies using eeganalparams
waveletFreqs = eeganalparams('freqs');
waveletWidth = eeganalparams('width');

%-- pre-allocate memory (to make sure it can be done!)
if PROCESS_CHANNELS_SEQUENTIALLY==0,  
    numChanPrealloc = numChannels;  
else
    numChanPrealloc = 1;
end
% #channels X #events X #freqs. X #timepoints = 4D array
arrayGB = numChanPrealloc * length(eventTrigger) * length(waveletFreqs) * DurationMS * 8 / 2^30; % 8 bytes per double, 2^30 bytes/GB
% initialize power matrices to make sure they can be stored
powerMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
powerMatZ = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
phaseMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);

clear defaultEEGfile subjDir eventEEGpath
% print statements for debugging and process checking
fprintf('\n');
fprintf('STEP 4 -- %d events to process for %s : %s', length(eventsTrig), subj, session);
fprintf('\n');
fprintf('The amount of RAM (GB) needed is: %d', arrayGB);
fprintf('\n\n');
disp('Variables to use here are:')
disp('powerMat, powerMatZ, phaseMat, numChanPrealloc')
disp('waveletFreqs, waveletWidth, ..')

disp(['size of matrices made are: '])
disp(size(powerMat))
fprintf('Number of preallocated channels are: %d', numChanPrealloc)
fprintf('\n');
fprintf('Duration of analysis is: %d', DurationMS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 5: Loop through the channels: extract, filter, processes, and save...   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Loop through and get all channels corresponding, and average
% loop through channels
% chanList = 46:96';
% chanStr = chanStr(46:96);
for iChan=1:numChannels
    powerMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    powerMatZ = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    phaseMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    
    
    thisChan = chanList(iChan,:);   % the channel to use in this loop (e.g. 48)
    thisChanStr = chanStr{iChan};
    strStart    = sprintf('\n STEP 5.%d -- Grab %d/%d: %s', iChan, iChan, numChannels, thisChanStr );  strStart(end+1:35)=' '; %buffer length so everything lines up
    fprintf('%s', strStart);       tic;

    %%- gete_ms: get the eegWaveV
    % eegwaveform for each event over the duration of time for a certain channel
    eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);

    % notch filter to eliminate 60 Hz noise
    fprintf(' [%.1f sec] --> notch filt', toc); tic;
    eegWaveV = buttfilt(eegWaveV,[59.5 60.5],resampledrate,'stop',1); %-filter is overkill: order 1 --> 25 dB drop (removing 5-15dB peak)
    
%     fs = 1000;
%     % robust spect parameters
%     alpha = 10;
%     window = 82;
%     
%     [xEst,freq,tWin,iter] = specPursuit(eegWaveV(1,:),fs,window,alpha);
% %     freq(1) = 0.01;
%     imagesc(tWin, log10(freq), 20*log10(abs(xEst)))
%     hold on; colormap(jet)
%     hCbar = colorbar('east');
%     set(hCbar,'ycolor',[1 1 1]*.1, 'fontsize', figFontAx-3, 'YAxisLocation', 'right')
%     set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
%     set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom

    %%- multiphasevec3: get the phase and power
    % power, phase matrices for events x frequency x duration of time for each channel
    fprintf(' [%.1f sec] --> freq decomp', toc); tic;
    [rawPhase,rawPow] = multiphasevec3(waveletFreqs,eegWaveV,resampledrate,waveletWidth);
    fprintf(' [%.1f sec] --> save', toc);  tic;
    fprintf('\n');
    
    % remove the leading/trailing buffer from the events we're interested in
    rawPow   = rawPow(:,:,BufferMS+1:end-BufferMS);
    rawPhase = rawPhase(:,:,BufferMS+1:end-BufferMS);
    eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
    eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
    if length(eegWaveT)<size(eegWaveV,2), % error check on time vs. voltage length
        fprintf('wave time length off'); 
        eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
    end

    % temp indicies
    iEv = 1:length(eventTrigger); % # of events
    iT  = 1:size(eegWaveV,2); % # of time points
    iF  = 1:length(waveletFreqs); % # of freqs.
    iChanSave = 1;

    % chan X event X freq X time
    % make power 10*log(power)
    powerMat(iChanSave,iEv,iF,iT) = 10*log10(rawPow);
    phaseMat(iChanSave,iEv,iF,iT) = rawPhase;

    % for each eegfile stem, z-score each channel and frequency
    fprintf(' [%.1f sec] --> z-score', toc);  tic;
    stemList = unique({eventTrigger.eegfile});
    for iStem=1:length(stemList),
        fprintf('.');
        iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));
        for iF = 1:length(waveletFreqs),
            allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
            mu = mean(allVal); stdev = std(allVal);
            
            % create the power matrix
            powerMatZ(iChanSave,iEvStem,iF,iT) = (powerMat(iChanSave,iEvStem,iF,iT)-mu)/stdev;
            
            if sum(isnan(powerMatZ(iChanSave,iEvStem,iF,iT)))>0
                keyboard;
            end
        end
    end
    fprintf(' [%.1f sec]', toc); tic;
    
    clear rawPow rawPhase
    disp('powerMatZ, powerMat and phaseMat are created')
    
    % normalized to 1 so all can be shifted to fit on a single plot
    eegWaveMeanSub  = eegWaveV-mean(mean(eegWaveV));   %double mean and double max because multiple events from same channel should be normalized together
    eegWaveShift    = (iChanSave-1)*2 + eegWaveMeanSub./max(max(abs(eegWaveMeanSub)));
    eegInstPow      = abs(eegWaveMeanSub).^2;
    eegInstPowShift = (iChanSave-1)*2 + eegInstPow./max(max(abs(eegInstPow)));
    wavesSft(iChanSave,iEv,iT) = eegWaveShift;

    % x-axis of time series
    waveT = eegWaveT;
    
    %% Save Data As Frequency/ProbeToVocalization Binned
    dataDir = 'condensed_data/';
    buffer_powerMatZ = squeeze(powerMatZ);
    
    %%- Call function to bin on freq. based on wavelet freqs. we have
    rangeFreqs = [freqBandAr.rangeF];
    rangeFreqs = reshape(rangeFreqs, 2, 7)';
    newPowerMatZ = freqBinSpectrogram(buffer_powerMatZ, rangeFreqs, waveletFreqs);
    
    buffer_powerMatZ = zeros(length(events), length(rangeFreqs));
    %%- Average From probe on -> vocalization [0: responseTime];
    responseTime = floor([events.responseTime]); % extract response Time's
    
    for i=1:length(events)
        buffer_powerMatZ(i,:) = mean(newPowerMatZ(i,:,2500:2500+responseTime(i)+1),3);
    end
    newPowerMatZ = buffer_powerMatZ
    
    %%- plotting for each trigger
    uniqueTrigType = unique(trigType);          % get all unique triggers (e.g. all probe words)
    numUniqueTrig  = length(uniqueTrigType);    % get length of unique triggers
    
    if length(chanList) == 1
        channel_num = chanList;
    else,
        channel_num = chanList(iChan);
    end
    
    clear buffer_powerMatZ 
    %%- Save this new power matrix Z
    data.uniqueTrigType = uniqueTrigType; % store all the unique trigger types
    data.trigType = trigType;             % store the trigger type per event
    data.powerMatZ = newPowerMatZ;        % save the condensed power Mat
    data.chanNum = thisChan;           % store the corresponding channel number
    data.chanStr = thisChanStr;               % the string name of the channel
    data.freqBandYtick = freqBandYticks;
    data.freqBandYlabel = freqBandYtickLabels;
    
    filename = strcat(dataDir, num2str(thisChan), '_', thisChanStr); 
    save(filename, 'data');
    
    clear data
    %% Finished looping through a certain Channel -> Save Data As Time/Frequency Binned
    %%- Here is where I save the data
    %%- i) condense data in freq/time domain and ii) save it
    dataDir = 'condensed_data/';
    if ~exist(dataDir), mkdir(dataDir); end
    
    %%- average data within time and within frequency
    
    % squeeze out the channel dimension of powerMatZ and time bin
    buffer_powerMatZ = squeeze(powerMatZ);
    %%- TIME BIN Before saving
    WinLength = 100; % 100 ms
    Overlap = 50;
    NumWins = size(buffer_powerMatZ,3) / (WinLength-Overlap) - 1;
    
    %%- Call function to bin on time based on winLength and Overlap and NumWins
    newPowerMatZ = timeBinSpectrogram(buffer_powerMatZ, NumWins, WinLength, Overlap);

    %%- Call function to bin on freq. based on 
    rangeFreqs = [freqBandAr.rangeF];
    rangeFreqs = reshape(rangeFreqs, 2, 7)';
    newPowerMatZ = freqBinSpectrogram(newPowerMatZ, rangeFreqs, waveletFreqs);
    
    timeZero = find(waveT==0,1)/Overlap + 1; % index of timezero in bins
    
    % create time vector that is binned and still centered at 0
    binnedWaveT = 1:size(newPowerMatZ,3) - timeZero;
    
    %%- plotting for each trigger
    uniqueTrigType = unique(trigType);          % get all unique triggers (e.g. all probe words)
    numUniqueTrig  = length(uniqueTrigType);    % get length of unique triggers
    
    if length(chanList) == 1
        channel_num = chanList;
    else,
        channel_num = chanList(iChan);
    end
    
    %% Create power mat data struct
    %save data that can plot evoked and spectrogram
    data.uniqueTrigType = uniqueTrigType; % store all the unique trigger types
    data.trigType = trigType;             % store the trigger type per event
    data.powerMatZ = newPowerMatZ;        % save the condensed power Mat
    data.waveT = waveT;             % save the binned Wave T
    data.waveletFreqs = waveletFreqs;     % wavelet frequencies used
    data.chanNum = thisChan;           % store the corresponding channel number
    data.chanStr = thisChanStr;               % the string name of the channel
    data.freqBandYtick = freqBandYticks;
    data.freqBandYlabel = freqBandYtickLabels;
    data.timeZero = timeZero;
    %%%%% ADD INDICES OF EVENTS?
    
    filename = strcat(dataDir, num2str(thisChan), '_', thisChanStr); 
    save(filename, 'data');
    
    %%- Clear useless vars
    clear buffer_powerMatZ WinLength Ovelap NumWins rangeFreqs data phaseMat powerMat powerMatZ
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 5i: Setup Vars to Plot  -------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print statements for debugging and process checking
fprintf('\n');
fprintf('STEP 5i -- %d events to plot evoked for ', length(eventsTrig));
fprintf('\n\n');

numChan = size(powerMatZ,1);
if numChan==1 && length(chanStr)~=1, %- numChan~=length(chanStr) means processed sequentially... make sure title refects correct channel
    %keyboard
    chanStr{1} = thisChanStr;
end

if length(chanList) == 1
    channel_num = chanList;
end

do_save = 0;
if do_save
    %%- TIME BIN Before saving
    WinLength = 100; % 100 ms
    Overlap = 50;
    NumWins = size(powerMatZ,4) / (WinLength-Overlap) - 1;

    timeZero = find(waveT==0,1)/Overlap + 1; % index of timezero in bins

    % initialize new power matrix
    newPowerMatZ = zeros(size(powerMatZ,1),size(powerMatZ,2),NumWins); 

    % window by time
    for i=1:size(powerMatZ,1),
        clear avgpowermat
        % loop through and bin the time domain into windows
        for j=1:NumWins %1:49:size(powerMatZ,3)-1
            % index through the time domain of the Z power matrix
            index = 1+(j-1)*50;

            if j == NumWins
                eventpowerMat = powerMatZ(i,:,index:end);
            else
                eventpowerMat = powerMatZ(i,:,index:index+100);
            end

            if ~exist('avgpowermat')
                avgpowermat = mean(eventpowerMat,3);
            else
                avgpowermat = [avgpowermat; mean(eventpowerMat,3)];
            end
        end
        % append avgpowermat to new matrix
        newPowerMatZ(i,:,:) = avgpowermat';
    end

    % create time vector that is binned and still centered at 0
    binnedWaveT = 1:size(newPowerMatZ,3) - timeZero;
    
    data.powerMatZ = newPowerMatZ;
    data.waveT = binnedWaveT;
    
    %%- Create power mat data struct
    %save data that can plot evoked and spectrogram
    data.uniqueTrigType = uniqueTrigType; % store all the unique trigger types
    data.trigType = trigType;             % store the trigger type per event
%     data.wavesSft = wavesSft;             % store teh eeg wave forms
%     data.waveT = waveT;                   % time points for the wave form
%     data.powerMatZ = powerMatZ;           % z-scored power
    data.waveletFreqs = waveletFreqs;     % wavelet frequencies used
    data.chanNum = channel_num;           % store the corresponding channel number
    data.chanStr = chanStr;               % the string name of the channel
    data.freqBandYtick = freqBandYticks;
    data.freqBandYlabel = freqBandYtickLabels;

    %%- Create raw data struct
    rawdata.waveletFreqs = waveletFreqs;    % store the wavelet freqs.
    rawdata.eegWaveV = eegWaveV;            % store the voltage time series per event
    rawdata.resampledrate = resampledrate;  % store resampling rate
    rawdata.waveletWidth = waveletWidth;    % width of wavelets
    rawdata.uniqueTrigType = uniqueTrigType;% the unique trigger types (e.g. brick, clock)
    rawdata.trigType = trigType;            % trigger per event
    rawdata.waveT = waveT;                  % the time points per eeg voltage
    rawdata.chanNum = channel_num;
    rawdata.chanStr = chanStr;
    rawdatafreqBandYtick = freqBandYticks;
    rawdata.freqBandYlabel = freqBandYtickLabels;

    chanString = chanStr{1};
    
    % call function to save raw data and data
    success = saveChannelData_powerAndRaw(data, rawdata, channel_num, chanString)
    clear data rawdata;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 6: Plot Evoked and Spectrogram For each trigger-------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numUniqueTrig>1,
    chanNum = numChan; % set the cahnnel number
    
    %-- Select Figure and open/hide --%
    thisFigNum = 1000+chanNum+FIG_OFFSET;
    if ~ishghandle(thisFigNum), 
        figure(thisFigNum);
    else
        set(0,'CurrentFigure', thisFigNum); 
    end
    if HIDE_FIGURES==1,     
        set(thisFigNum,'visible','off');    
    end 
    clf
    set(gcf,'color','w')
    cTrigAr = 'rkbgmcy';
    
    %-- figure title
    titleStr = sprintf('%s : %s : channel %s number(%s)', subj, session, chanStr{chanNum}, num2str(chanList));
    titleStr(find(titleStr=='_'))=' ';
    title(titleStr, 'fontsize',20)
    
    for thisTrig = 1:numUniqueTrig,
        % find all event indices for thisTrig (e.g. brickprobe)
        iTrig = find(strcmp(trigType,uniqueTrigType(thisTrig)));
        chanNum = iChan; % set the current channel number
        
        %%----------------- plot event types--------------------------------------------------------
        subplot(2+numUniqueTrig,1,1)
        axTask = gca;
            
        % make a struct to hold all the events for this trigger
        eventTypes   = eventsTrig(iTrig);
        yValues  = thisTrig; % increment sequentially... don't use value of trigType, which can jump
        
        % ? depends on mstimeEnd which there is no field inside events
%         tOffset  = trigZeroMS(iTrig);
              
        %-- event start and stop time
        % ? no .msDuration though either?
%         eTimeOnS  = ([events.mstime]-tOffset)/1000;
%         eTimeOffS = ([events.mstime]+[events.msDuration]-tOffset)/1000;
%         
%         % all events get a blue + to indicate start
%         hP = plot(eTimeOnS,yValues,'b.','MarkerSize',10); hold on;  
%         % all events get a line indicating duration
%         hL = line(median([eTimeOnS; eTimeOffS],2), [yValues; yValues],'Color',cTrig,'LineWidth',4);
%         
%         %-- next (meta) event start time
%         % ? no .msMetaEventNext
%         eTimeEnxtS = ([events.msMetaEventNext]-tOffset)/1000;   % time that next meta event starts
%         hP = plot(eTimeEnxtS,yValues-.1,'k>','MarkerSize',6);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%------------------ STEP 6i: Plot Evoked -------------------------%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % print statements for debugging and process checking
        fprintf('Looking at %s probe word', uniqueTrigType{thisTrig});
        fprintf('\n');
        fprintf('STEP 6i -- %d events to plot evoked for %s', length(iTrig), chanStr{chanNum});
        fprintf('\n\n');
    
        %-compute the evoked response for each trigger type
        wavePlot_allEvents = wavesSft(:,iTrig,:); % get waveplot for the specific triggers
        clear evokedRespMu evokedRespSEM;
        % reduce dimensionality for plotting (average across all events)
        evokedRespMu(1:size(wavePlot_allEvents,3))  = mean(wavePlot_allEvents(chanNum,:,:),2);
        evokedRespSEM(1:size(wavePlot_allEvents,3)) = std(wavePlot_allEvents(chanNum,:,:),0,2)./sqrt(size(wavePlot_allEvents,2));
    
        cTrigAr = 'rkbgmcy';
        cTrig = cTrigAr(thisTrig); % set color scheme
        
        % actually plot evoked mean
        subplot(2+numUniqueTrig,1,2) % put all evoked responses on 1 plot
        axWave = gca;
        hEvoked = plot(waveT, evokedRespMu, 'k-'); hold on
        set(hEvoked,'Color',cTrig,'LineWidth',2)
        axis tight;
        waveYlabel = 'Evoked Response (norm)';
        ylabel(waveYlabel,'fontsize',figFontAx)
        
%         ii    = length(evokedRespMu):-1:1;
%         fillX = [waveT waveT(ii)];
%         fillY = [evokedRespMu+evokedRespSEM evokedRespMu(ii)-evokedRespSEM(ii)];
%         hFill = fill(fillX,fillY,cTrig);
%         hLegIn(thisTrig)  = hFill;
%         hLegStr{thisTrig} = sprintf('%s [%d ev] .',uniqueTrigType(thisTrig),size(wavePlot_allEvents,2));  %
        
        clear wavePlot_allEvents;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%------------------ STEP 6ii: Plot Spectrogram -------------------------%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % print statements for debugging and process checking
        fprintf('\n');
        fprintf('STEP 6ii -- %d events to plot evoked for %s', length(iTrig), chanStr{chanNum});
        fprintf('\n\n');
        
        % get the power for corresponding events
        thisPowMat = powerMatZ(chanNum,iTrig,:,:); 
        if thisTrig==1, 
            cAx = [-.8 .8]; 
            if length(iTrig)<4,
                cAx=cAx*3; 
            elseif length(iTrig)>100, 
                cAx=cAx/2; 
            end;   
        end;   % z-scored power... set cAx for first panel and keep same for all others

        % reduce the dimensionality of the matrix for plotting - average
        % across all events
        powPlot = mean(thisPowMat,2); 
        titleStr = sprintf('mean power: chan %s, %d events', chanStr{chanNum}, size(thisPowMat,2));
        powPlot = squeeze(powPlot); % squeeze out the singleton dimension
        
         % actually plot the spectrogram
        subplot(2+numUniqueTrig,1,2+thisTrig)
        axSpec(thisTrig)=gca;
        hImg    = imagesc(waveT,log10(waveletFreqs),powPlot); 
        hold on;  colormap(jet);

        hCbar = colorbar('east');
        set(hCbar,'ycolor',[1 1 1]*.1, 'fontsize', figFontAx-3, 'YAxisLocation', 'right')
            
        % set the heat map settings
        set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
        set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
        set(gca,'fontsize',figFontAx+3)
        set(gca,'XTick',[],'Box','off');
        title(uniqueTrigType(thisTrig), 'fontsize',20)
%         set(gca,'clim',[-1 1]*0.5)      
    end % end of loop through unique triggers  
end

% set legend at end
subplot(1+numUniqueTrig,1,1) % put all evoked responses on 1 plot
axWave = gca;
legend(uniqueTrigType);