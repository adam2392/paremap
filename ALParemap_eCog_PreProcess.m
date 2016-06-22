%         -- preprocess data with notch filter, wavelet transform
%         -- select subset of behavioral events (filter out incorrect
%         responses)
%         -- get waveform, power, phase for a time window around the event of interest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subject, timeLock


function ALParemap_eCog_PreProcess(subj, typeTransform, timeLock, referenceType, winSize, stepSize)
% clear all;
close all;
clc;

% subj = 'NIH034';
% timeLock = 'vocalization';
% referenceType = 'global';
% winSize = 500;
% stepSize = 100;
% typeTransform = 'multitaper';

USE_CHAN_SUBSET = 0; % 0=all channels, 1=process the subset
%% PARAMETERS FOR RUNNING PREPROCESS
expected_timeLocks = {'vocalization', 'matchword', 'probeword'};
expected_transforms = {'morlet', 'multitaper'};
REF_TYPES = {'noreref', 'bipolar', 'global'};

DEBUG = 1;
if ~ismember(timeLock, expected_timeLocks)
    disp('timeLock should be vocalization, matchword, or probeword');
end
if ~ismember(referenceType, REF_TYPES)
    disp('reference types are noreref, bipolar, or global');
end
if ~ismember(typeTransform, expected_transforms)
    disp('transform types are morlet, multitaper');
end

% REFERENCE ELECTRODE
THIS_REF_TYPE = referenceType; 
% FILTERING OPTIONS
BP_FILTER_RAW                 = 1;  %-0 or 1: apply a bandpass filter to the raw traces (1-499 hz)

% array of frequency bands
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

%% LOAD EVENTS STRUCT AND SET DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
eegRootDirJhu = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

% Either go through all the sessions, or a specific session
disp('STEP 1: Going through all sessions')
session = 'Meta Session [all]';
behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap');

subjDir = fullfileEEG(eegRootDir,subj); % directory to subject (e.g. NIH034)
docsDir = fullfileEEG(subjDir,'docs');  % directory to the docs (electordes.m, tagNames.txt, etc.)
talDir  = fullfileEEG(subjDir,'tal');
defaultEEGfile = fullfileEEG('/Volumes/Shares/FRNU/data/eeg/',subj,'/eeg.reref/');  % default event eegfile fields point here... switch to local before loading

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

%%- GET CORRECT EVENTS ONLY
% POST MODIFY EVENTS based on fields we want (e.g. is it correct or not)?
correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);

%% LOAD CHANNEL LIST AND INFORMATION
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
        fprintf('Bipolar referencing');
        chanFile      = [talDir '/leads_bp.txt'];
        [chan1 chan2] = textread(chanFile,'%d%*c%d');
        chanList      = [chan1 chan2];
        for iChan=1:size(chanList,1),
            %    chanStr{iChan} = sprintf('%d-%d (%s-%s)', chan1(iChan), chan2(iChan), chanTags{find(chanNums==chan1(iChan))}, chanTags{find(chanNums==chan2(iChan))} );
            chanStr{iChan} = sprintf('%s-%s', chanTags{find(chanNums==chan1(iChan))}, chanTags{find(chanNums==chan2(iChan))} );
        end
        chanRefs      = [];
        eventEEGpath  = '/eeg.reref/';
    case 'global' % look at global electrodes / monopolar
        fprintf('STEP 1: Using Global referencing\n');
        chanFile      = [talDir '/leads.txt'];
        chanList      = textread(chanFile,'%d'); % read in the list of channels nums

        % set the names for each channel
        for iChan=1:size(chanList,1),
            chanStr{iChan} = sprintf('%s-global', chanTags{find(chanNums==chanList(iChan))} );
        end
        eventEEGpath  = '/eeg.reref/';
    otherwise
        fprintf('Error, no referencing scheme selected');
end

iChanListSub  = 33:96;            %G1, G2, LF1, AST1,
%%- select all channels, or part of the subset of channels
if USE_CHAN_SUBSET==0,
    iChanList = 1:size(chanList,1);  %all possible channels
else
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

%% FILTER AND PREPROCESS PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 3: DATA INPUT TO GETE_MS, MULTIPHASEVEC3 AND ZSCORE
%%------------------      AND SET UP POWER, POWERZ AND PHASE MATRICS        ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- Input to gete_ms
%%- Dependent only on eventsTriggerXlim: These stay the same regardless of how we process events
eventTrigger = events;

% offset to synchronize with vocalization
if strcmp(timeLock, 'vocalization'),
    for iEvent=1:length(eventTrigger),
        eventTrigger(iEvent).mstime = eventTrigger(iEvent).mstime + eventTrigger(iEvent).responseTime;
        eventTrigger(iEvent).eegoffset = eventTrigger(iEvent).eegoffset + round(eventTrigger(iEvent).responseTime);
    end
    LOWERTIME = -2;
    UPPERTIME = 3;
    disp('VOCALIZATION');
elseif strcmp(timeLock, 'matchword'),
    for iEvent=1:length(eventTrigger),
        eventTrigger(iEvent).eegoffset = eventTrigger(iEvent).eegoffset + round(eventTrigger(iEvent).matchOnTime - eventTrigger(iEvent).mstime);
        eventTrigger(iEvent).mstime = eventTrigger(iEvent).mstime + (eventTrigger(iEvent).matchOnTime - eventTrigger(iEvent).mstime);
    end
    LOWERTIME = -4;
    UPPERTIME = 1;
    disp('MATCHWORD');
elseif strcmp(timeLock, 'probeword')
    % Settings for probewordon synchronization
    LOWERTIME = -1;
    UPPERTIME = 4;
    disp('PROBE WORD ON');
else
    error('not set correctly.');
end

eventsTriggerXlim = [LOWERTIME UPPERTIME]; % range of time to get data from (-2 seconds to 5 seconds after mstime (probeWordOn)) 
eventOffsetMS   = eventsTriggerXlim(1)*1000;      % positive = after event time; negative = before event time
eventDurationMS = diff(eventsTriggerXlim)*1000;   % duration includes offset (i.e., if offset -500 and duration 1000, only 500 ms post event will be prsented)

OffsetMS        = eventOffsetMS;     % positive = after event time; negative = before event time
DurationMS      = eventDurationMS;   % duration includes offset (i.e., if offset -500 and duration 1000, only 500 ms post event will be prsented)
BufferMS        = 1000;              % grab excess data before/after event window so filters don't have edge effect
resampledrate   = 1000;              % don't resample... keep the 1kHz native sample rate

%%- apply a bandpass filter raw data? (i.e. pre-filter the wave?)
if BP_FILTER_RAW==1,
    preFiltFreq      = [1 499];   %[1 499] [2 250]; first bandpass filter data from 1-499 Hz
    preFiltType      = 'bandpass';
    preFiltOrder     = 2;
    preFiltStr       = sprintf('%s filter raw; %.1f - %.1f Hz',preFiltType,preFiltFreq);
    preFiltStrShort  = '_BPfilt';
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
    eventTrigger(iEvent).eegfile = regexprep(eventTrigger(iEvent).eegfile, defaultEEGfile, fullfileEEG(subjDir,eventEEGpath));
end

%%- gets the range of frequencies using eeganalparams
waveletFreqs = eeganalparams('freqs');
waveletWidth = eeganalparams('width');

%-- pre-allocate memory (to make sure it can be done!)
% #channels X #events X #freqs. X #timepoints = 4D array
arrayGB = length(eventTrigger) * length(waveletFreqs) * DurationMS * 8 / 2^30; % 8 bytes per double, 2^30 bytes/GB
clear defaultEEGfile subjDir eventEEGpath
% print statements for debugging and process checking
fprintf('\n');
fprintf('STEP 4 -- %d events to process for %s : %s', length(eventTrigger), subj, session);
fprintf('\n');
fprintf('The amount of RAM (GB) needed is: %d', arrayGB);
fprintf('\n\n');
disp('Variables to use here are:')
disp('powerMat, powerMatZ, phaseMat')
disp('waveletFreqs, waveletWidth, ..')
fprintf('Duration of analysis is: %d\n', DurationMS)

%% EXTRACT, FILTER, PROCESS AND SAVE PREPROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 5: Loop through the channels: extract, filter, processes, and save...   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROBUST_SPEC = 0; % carry out robust spectrogram on data instead
ALPHA = 300;    % l1 regularization parameter
WINDOW = 200;   % window size for spectrotemporal pursuit
FS = 1000;      % sampling frequency of data

SAVE = 1;       % save data?
for iChan=1:numChannels
    powerMat  = zeros(length(eventTrigger), length(waveletFreqs), DurationMS);
    powerMatZ = zeros(length(eventTrigger), length(waveletFreqs), DurationMS);
    phaseMat  = zeros(length(eventTrigger), length(waveletFreqs), DurationMS);
    
    thisChan = chanList(iChan,:);   % the channel to use in this loop (e.g. 48)
    thisChanStr = chanStr{iChan};
    strStart    = sprintf('\n STEP 5.%d -- Grab %d/%d: %s', iChan, iChan, numChannels, thisChanStr );  strStart(end+1:35)=' '; %buffer length so everything lines up
    fprintf('%s', strStart);   tic;
    
    %% SETTINGS FOR SAVING PROCESSED DATA 
    chanFileName = strcat(num2str(thisChan), '_', thisChanStr);

    % data directories to save data into
    workDir = '/Users/liaj/Documents/MATLAB/paremap';
    homeDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/';
    jhuDir = '/home/adamli/paremap/';

    % Determine which directory we're working with automatically
    if     ~isempty(dir(workDir)), rootDir = workDir;
    elseif ~isempty(dir(homeDir)), rootDir = homeDir;
    elseif ~isempty(dir(jhuDir)), rootDir = jhuDir;
    else   error('Neither Work nor Home EEG directories exist! Exiting'); end

    subjDir = strcat('condensed_data_', subj);
    typeTransformDir = fullfile(rootDir, subjDir, strcat(typeTransform, '_', referenceType));
    responseDir = fullfile(typeTransformDir, timeLock);

    %% EXTRACT WAVEFORM, TRANSFORM, LOAD, SAVE
    %%- gete_ms: get the eegWaveV
    % eegwaveform for each event over the duration of time for a certain channel
    eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);

    % notch filter to eliminate 60 Hz noise
    fprintf(' [%.1f sec] --> notch filt\n', toc); tic;
    eegWaveV = buttfilt(eegWaveV,[59.5 60.5],resampledrate,'stop',1); %-filter is overkill: order 1 --> 25 dB drop (removing 5-15dB peak)
        
    if strcmp(typeTransform, 'morlet') % OPTION 1: perform wavelet spectral analysis
        %%- multiphasevec3: get the phase and power
        % power, phase matrices for events x frequency x duration of time for each channel
        fprintf(' [%.1f sec] --> freq decomp', toc); tic;
        [rawPhase,rawPow] = multiphasevec3(waveletFreqs,eegWaveV,resampledrate,waveletWidth);
        fprintf(' [%.1f sec] --> save', toc);  tic;
        fprintf('\n');

        %%- REMOVE LEADING/TRAILING buffer areas from power, phase,
        %%eegWave, timeVector
        rawPow   = rawPow(:,:,BufferMS+1:end-BufferMS);
        rawPhase = rawPhase(:,:,BufferMS+1:end-BufferMS);
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
        eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
        if length(eegWaveT)<size(eegWaveV,2), % error check on time vs. voltage length
            fprintf('wave time length off'); 
            eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
        end
        % x-axis of time series
        waveT = eegWaveT;

        % temp indicies
        iEv = 1:length(eventTrigger); % # of events
        iT  = 1:size(eegWaveV,2); % # of time points
        iF  = 1:length(waveletFreqs); % # of freqs.

        % chan X event X freq X time
        % make power 10*log(power)
        powerMat(iEv,iF,iT) = 10*log10(rawPow);
        phaseMat(iEv,iF,iT) = rawPhase;
        
        % set two paramters from robust spectrotemp to 0 
        tWin = 0;
        freq = 0;
        freqs = waveletFreqs;
        
        % create vector of the actual seconds in time axis for the powerMat
        % (since its time binned)...
        OVERLAP = stepSize;
        WINSIZE = winSize;
        if tWin == 0, % if not set yet
            tWin = (LOWERTIME) :OVERLAP/FS: (UPPERTIME)-WINSIZE/FS;
        end
        timeZero = abs(0-(LOWERTIME))/(OVERLAP/FS) - (WINSIZE - OVERLAP)/OVERLAP;

        
        fprintf(' [%.1f sec]', toc); tic;
        clear rawPow rawPhase
        disp('powerMatZ, powerMat and phaseMat are created')
    elseif strcmp(typeTransform, 'multitaper')
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
        eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
     
        Fs = 1000;
        T = winSize/1000;
        overlap = stepSize/1000;
        mtBandWidth = 4;
        mtFreqs = [];

        %%- multitaper FFT 
        [rawPowBase, freqs_FFT, t_sec,rawPhaseBase] = eeg_mtwelch2(eegWaveV, Fs, T, overlap, mtBandWidth, mtFreqs, 'eigen');
        fprintf(' [multitaper%.1f|',toc); tic;
        freq = freqs_FFT;
        
        powerMat = 10*log10(rawPowBase);
       
        freqs = freqs_FFT;
        
        % temp indicies
        iEv = 1:length(eventTrigger); % # of events
        iT  = 1:size(powerMat,3); % # of time points
        iF  = 1:length(freqs); % # of freqs.
        
        tWin = t_sec + LOWERTIME;
        [~, timeZero] = min(abs(tWin(:,2) + tWin(:,1)));
%         rawPSD = zeros(size(rawPowBase));
%         rawPhase = zeros(size(rawPhaseBase));
%         for iEvent=1:size(eegWaveV, 1),
%             X = eegWaveV(iEvent, :);
%             [~, PS, PSD, freqs_FFT, all_phase, t_sec] = fft_spectral(X,Fs,T,overlap, 'multitaper', 4, 'weights', 'eigen');
%             rawPSD(iEvent,:,:) = PSD;
%             
%             if ~isempty(find(squeeze(rawPowBase(iEvent,:,:)) == PSD))
%                 disp('wtf')
%             end
%         end
%         fprintf(' [multitaper %.1f|',toc); tic;
    elseif strcmp(typeTransform, 'robust_spec') % OPTION 2: Perform robust spectrotemporal pursuit instead 
        %%- REMOVE LEADING/TRAILING BUFFER REGION FOR EEGWAVE
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
        eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
        if length(eegWaveT)<size(eegWaveV,2), % error check on time vs. voltage length
            fprintf('wave time length off'); 
            eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
        end
        
        %%- ROBUST SPECTROTERMPORAL PURSUIT PARAMETERS
        % initialize arrays for storing 
        powerMat = zeros(length(eventTrigger), WINDOW/2, length(eegWaveV)/WINDOW); % #EVENTS X #FREQ X #TIME
        powerMatZ = zeros(length(eventTrigger), WINDOW/2, length(eegWaveV)/WINDOW);
        
        % generate robust spectrogram for each event
        tic;
        for iEv=1:length(eventTrigger),
            [xEst,freq,tWin,iter] = specPursuit(eegWaveV(iEv,:), FS, WINDOW, ALPHA);
            size(xEst)
            xEst = 20*log10(abs(xEst)); % convert to power spectrum
            xEst = reshape(xEst, 1, size(xEst,1), size(xEst,2));
            powerMat(iEv,:,:) = xEst;
        end
        fprintf(' [%.1f sec] --> robust spect pursuit', toc);
        
        % reset tWin var to reflect over seconds
        tWin = tWin - timeZero*1000;
        
        %%- Zscore all power for robust spectral pursuit
        fprintf(' [%.1f sec] --> z-score robust spec pursuit', toc); tic;
        for iEvent =1:size(powerMat,1),
            for iF = 1:size(powerMat,2),
                fixedVal = powerMat(iEvent, iF, 1:4); %allVal for particular chan and freq
                mu = mean(fixedVal); stdev = std(fixedVal);

                % create the power matrix
                powerMatZ(iEvent, iF, :) = (powerMat(iEvent, iF, :)-mu)/stdev;

                if sum(isnan(powerMatZ(iEvent,iF,:)))>0
                    keyboard;
                end
            end
        end
        fprintf(' [%.1f sec]', toc);
    end
    
    %     for each eegfile stem, z-score each channel and frequency
    fprintf(' [%.1f sec] --> z-score', toc);  tic;
    stemList = unique({eventTrigger.eegfile});
    powerMatZ = zeros(size(powerMat));
    
    % indices of the powerMat to Z-score wrt
    for iStem=1:length(stemList),
        fprintf('.');
        iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));

%             disp(['Length of session eeg file: ', length(iEvStem)])
        for iF = 1:length(freqs),
            allVal = reshape(squeeze(powerMat(iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
            mu = mean(allVal); stdev = std(allVal);

            % create the power matrix
            powerMatZ(iEvStem,iF,iT) = (powerMat(iEvStem,iF,iT)-mu)/stdev;

            if sum(isnan(powerMatZ(iEvStem,iF,iT)))>0
                keyboard;
            end
        end
    end
    
    if strcmp(typeTransform, 'morlet')
        %% TIME BIN POWERMATZ WITH WINDOWSIZE AND OVERLAP
        addpath('./m_oldAnalysis_anovaANDsinglechannel/');
        WINDOWSIZE = 500; % in milliseconds
        OVERLAP = 100;    % in milliseconds
        powerMatZ = timeBinSpectrogram(powerMatZ, WINDOWSIZE, OVERLAP);

        if DEBUG,
            size(powerMatZ)
        end

        %% FREQUENCY BIN WITH FREQUENCY BANDS
        rangeFreqs = reshape([freqBandAr.rangeF], 2, 7)';
        waveletFreqs = waveletFreqs;
        powerMatZ = freqBinSpectrogram(powerMatZ, rangeFreqs, waveletFreqs);

        if DEBUG,
            size(powerMatZ)
        end
    elseif strcmp(typeTransform, 'multitaper')
        %%- FREQUENCY BIN
    end
    
    clear powerMat
    if DEBUG,
        size(powerMatZ)
    end
    
    %%- SAVE 01: ENTIRE DATASET PER CHANNEL
    targetWords = unique({events.targetWord});
    for iTarget=1:length(targetWords)
        THIS_TARGET = targetWords{iTarget};
        
        eventIndices = find(strcmp({events.targetWord}, THIS_TARGET));
        thisPowMat = powerMatZ(eventIndices,:,:);
        size(thisPowMat)
        
        if SAVE,
            %%- Save this new power matrix Z-scored into data .mat file
            data.targetWords = THIS_TARGET;                 % the target words for all events in this struct
            data.eegWaveV = eegWaveV(eventIndices,:);       % eeg wave form
            data.eegWaveT = eegWaveT;                       % time series for eeg voltage wave
            data.chanNum = thisChan;                        % store the corresponding channel number
            data.chanStr = thisChanStr;                     % the string name of the channel
            data.freqBandYtick = 1:length(freqBandYticks);            % store frequency bands if using wavelet transform
            data.freqBandYlabel = {freqBandAr.name};
            data.descriptor = 'Initial processing -2 seconds to 4 seconds after VOCALIZATION. Time binned with 500ms window and 100ms overlap';
            data.timeZero = timeZero; %ceil((TIMEZERO-LOWERTIME)/OVERLAP);
            data.powerMatZ = thisPowMat;            % save the condensed power Mat Z-scored
%             data.waveT = tWin;                      % ROBUSTSPECT: save the binned Wave T
%             data.freq = freq;                       % ROBUSTSPECT: save the frequency points

            fileDir = strcat(responseDir, '_targetWords/', THIS_TARGET);
            chanFilePath = fullfile(fileDir, chanFileName);
            
            if ~exist(fileDir)
                mkdir(fileDir);
            end
            save(chanFilePath, 'data');            
        end
    end % loop through targetwords

    %%- SAVE 02: VOCALIZED WORDS PER SESSION/BLOCK
    % SPLIT INTO SESSIONS AND BLOCKS
    subjSessions = unique({events.sessionName}); % e.g. sessions 0, 1, 2
    subjBlocks = unique({events.blocknumber});   % e.g. blocks 0,1,2,3,4,5
    
    for iSesh=1:length(subjSessions),
        for iBlock=1:length(subjBlocks),
            sessionBlockIndices = strcmp({events.sessionName}, subjSessions(iSesh)) & ...
                                    strcmp({events.blocknumber}, subjBlocks(iBlock));
             
            %%- ONLY ANALYZE TARGET WORDS -> VOCALIZATION OF 'S' SOUNDING
%             targetWords = unique({events(sessionBlockIndices).targetWord});
%             
%             for iTarget=1:length(targetWords)
%                 THIS_TARGET = targetWords{iTarget};
%                 eventIndices = find(strcmp({events.targetWord}, THIS_TARGET) & ...
%                                     sessionBlockIndices);
%                 
%                 sessionBlockVocalWordEvents = events(eventIndices); 
%                 blockNum = unique({sessionBlockVocalWordEvents.blocknumber});
%                 sessionNumber = sessionBlockVocalWordEvents(1).sessionNum;
%                     
%                 thisPowMat = powerMatZ(eventIndices,:,:);
%                 
%                 %% SAVE PROCESSED DATA IN A MATLAB STRUCT
%                 if SAVE,
%                     %%- Save this new power matrix Z-scored into data .mat file
%                     data.targetWords = THIS_TARGET;                 % the target words for all events in this struct
%                     data.sessionNum = sessionNumber;                % the session number
%                     data.blockNum = blockNum;                       % the block number
%                     data.eegWaveV = eegWaveV(eventIndices,:);       % eeg wave form
%                     data.eegWaveT = eegWaveT;                       % time series for eeg voltage wave
%                     data.chanNum = thisChan;                        % store the corresponding channel number
%                     data.chanStr = thisChanStr;                     % the string name of the channel
%                     data.freqBandYtick = 1:length(freqBandYticks);            % store frequency bands if using wavelet transform
%                     data.freqBandYlabel = {freqBandAr.name};
%                     data.descriptor = 'Initial processing -2 seconds to 4 seconds after VOCALIZATION. Time binned with 500ms window and 100ms overlap';
%                     data.timeZero = timeZero; %ceil((TIMEZERO-LOWERTIME)/OVERLAP);
%                     data.vocalization = data.timeZero + ceil([sessionBlockVocalWordEvents.responseTime]/OVERLAP);
%                     data.powerMatZ = thisPowMat;            % save the condensed power Mat Z-scored
%                     data.waveT = tWin;                      % ROBUSTSPECT: save the binned Wave T
%                     data.freq = freq;                       % ROBUSTSPECT: save the frequency points
%                     
%                     %%- SAVING DIR PARAMETERS
%                     if ROBUST_SPEC,
%                         TYPE_SPECT = 'robust_spec';
%                     else
%                         TYPE_SPECT = 'morlet_spec';
%                     end
%                     if VOCALIZATION
%                         TYPE_SPECT = strcat(TYPE_SPECT, '_vocalization');
%                     elseif MATCHWORD
%                         TYPE_SPECT = strcat(TYPE_SPECT, '_matchword');
%                     end
% 
%                     chanFileName = strcat(num2str(thisChan), '_', thisChanStr, '_', TYPE_SPECT);
% 
%                     % data directories to save data into
%                     workDir = '/Users/liaj/Documents/MATLAB/paremap';
%                     homeDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/';
%                     jhuDir = '/home/adamli/paremap/';
% 
%                     % Determine which directory we're working with automatically
%                     if     length(dir(workDir))>0, rootDir = workDir;
%                     elseif length(dir(homeDir))>0, rootDir = homeDir;
%                     elseif length(dir(jhuDir))>0, rootDir = jhuDir;
%                     else   error('Neither Work nor Home EEG directories exist! Exiting'); end
% 
%                     dataDir = strcat('condensed_data_', subj);
%                     typeTransformDir = fullfile(rootDir, dataDir, TYPE_SPECT);
%                     fileDir = fullfile(typeTransformDir, subjSessions{iSesh}, subjBlocks{iBlock}, THIS_TARGET);
%                     chanFilePath = fullfile(fileDir, chanFileName);; 
% 
%                     if ~exist(fileDir)
%                         mkdir(fileDir);
%                     end
%                     save(chanFilePath, 'data');            
%                 end
%             end
            
            %%- SAVE 03: WORD PAIRS PER SESSION/BLOCK
            %%- CARRY OUT REGULAR ANALYSIS ON WORD PAIR GROUPS
            %%- LOOP THROUGH PROBEWORDS
            probeWords = unique({events(sessionBlockIndices).probeWord});
            for iProbe=1:length(probeWords),
                THIS_PROBE = probeWords{iProbe};
                
                % events for this probeWord and their targetWords
                probeIndices = strcmp({events.probeWord}, THIS_PROBE);
                tempEvents = events(probeIndices & sessionBlockIndices);
                targetWords = unique({tempEvents.targetWord});
                
                %%- LOOP THROUGH TARGETWORDS FOR EACH PROBEWORD
                for iTarget=1:length(targetWords),
                    THIS_TARGET = targetWords{iTarget};
                    
                    % match probe, target, session and block
                    eventIndices = find(strcmp({events.probeWord}, THIS_PROBE) & ...
                                    strcmp({events.targetWord}, THIS_TARGET) & ...
                                    sessionBlockIndices);
                    sessionBlockWordPairEvents = events(eventIndices);
                    
                    blockNum = unique({sessionBlockWordPairEvents.blocknumber});
                    sessionNumber = sessionBlockWordPairEvents(1).sessionNum;
                    
                    thisPowMat = powerMatZ(eventIndices,:,:);
                    
                    %% SAVE PROCESSED DATA IN A MATLAB STRUCT
                    if SAVE,
                        %%- Save this new power matrix Z-scored into data .mat file
                        data.probeWords = THIS_PROBE;                   % the probe words for all events in this struct
                        data.targetWords = THIS_TARGET;                 % the target words for all events in this struct
                        data.sessionNum = sessionNumber;                % the session number
                        data.blockNum = blockNum;                       % the block number
                        data.eegWaveV = eegWaveV(eventIndices,:);       % eeg wave form
                        data.eegWaveT = eegWaveT;                       % time series for eeg voltage wave
                        data.chanNum = thisChan;                        % store the corresponding channel number
                        data.chanStr = thisChanStr;                     % the string name of the channel
                        data.freqBandYtick = 1:length(freqBandYticks);            % store frequency bands if using wavelet transform
                        data.freqBandYlabel = {freqBandAr.name};
                        data.descriptor = 'Initial processing -2 seconds to 3 seconds after VOCALIZATION. Time binned with 500ms window and 100ms overlap';
                        data.timeZero = timeZero; %ceil((TIMEZERO-LOWERTIME)/OVERLAP);
                        data.powerMatZ = thisPowMat;            % save the condensed power Mat Z-scored
                        data.waveT = tWin;                      % ROBUSTSPECT: save the binned Wave T
                        data.freq = freq;                       % ROBUSTSPECT: save the frequency points

                        %%- SAVING DIR PARAMETERS
                        chanFileName = strcat(num2str(thisChan), '_', thisChanStr);
                        wordpair_name = strcat(THIS_PROBE, '_', THIS_TARGET);
                        fileDir = fullfile(responseDir, subjSessions{iSesh}, subjBlocks{iBlock}, wordpair_name);
                        chanFilePath = fullfile(fileDir, chanFileName);; 

                        if ~exist(fileDir)
                            mkdir(fileDir);
                        end
                        save(chanFilePath, 'data');            
                    end
                end % loop through target
            end % loop through probe
        end % loop through block 
    end % loop through session
end
end