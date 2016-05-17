%         -- preprocess data with notch filter, wavelet transform
%         -- select subset of behavioral events (filter out incorrect
%         responses)
%         -- get waveform, power, phase for a time window around the event of interest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% PARAMETERS FOR RUNNING PREPROCESS
subj = 'NIH034';
sessNum = [0, 1, 2];
DEBUG = 1;

REF_TYPES = {'noreref', 'bipolar', 'global'};
THIS_REF_TYPE = REF_TYPES{3}; 

USE_CHAN_SUBSET = 0; % 0=all channels, 1=process the subset

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

% FILTERING OPTIONS
BP_FILTER_RAW                 = 1;  %-0 or 1: apply a bandpass filter to the raw traces (1-499 hz)
PROCESS_CHANNELS_SEQUENTIALLY = 1;  %0 or 1:  0 means extract all at once, 1 means sequentially

%% LOAD EVENTS STRUCT AND SET DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
% eegRootDirHome = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

% Either go through all the sessions, or a specific session
if sessNum == -1 | length(sessNum)>1, % all sessions
    disp('STEP 1: Going through all sessions')
    session = 'Meta Session [all]';
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap');
    sessStr = '[all]';
else                                  % one session
    disp('STEP 1: Going through one session')
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap/', session);
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
    case 'global' % look at global electrodes
        fprintf('STEP 1: Using Global referencing\n');
        chanFile      = [talDir '/leads.txt'];
        chanList      = textread(chanFile,'%d'); % read in the list of channels nums
        
        % set the names for each channel
        for iChan=1:size(chanList,1),
            chanStr{iChan} = sprintf('%s-global', chanTags{find(chanNums==chanList(iChan))} );
        end
        eventEEGpath  = '/eeg.reref/';
        
        iChanListSub  = 2:96;            %G1, G2, LF1, AST1,
    otherwise
        fprintf('Error, no referencing scheme selected');
end

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
LOWERTIME = -1;
UPPERTIME = 5;
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
fprintf('STEP 4 -- %d events to process for %s : %s', length(eventTrigger), subj, session);
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
    powerMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    powerMatZ = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    phaseMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    
    thisChan = chanList(iChan,:);   % the channel to use in this loop (e.g. 48)
    thisChanStr = chanStr{iChan};
    strStart    = sprintf('\n STEP 5.%d -- Grab %d/%d: %s', iChan, iChan, numChannels, thisChanStr );  strStart(end+1:35)=' '; %buffer length so everything lines up
    fprintf('%s', strStart);   tic;
    
    %%- gete_ms: get the eegWaveV
    % eegwaveform for each event over the duration of time for a certain channel
    eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);

    % notch filter to eliminate 60 Hz noise
    fprintf(' [%.1f sec] --> notch filt\n', toc); tic;
    eegWaveV = buttfilt(eegWaveV,[59.5 60.5],resampledrate,'stop',1); %-filter is overkill: order 1 --> 25 dB drop (removing 5-15dB peak)
    
    if ~ROBUST_SPEC % OPTION 1: perform wavelet spectral analysis
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

        % temp indicies
        iEv = 1:length(eventTrigger); % # of events
        iT  = 1:size(eegWaveV,2); % # of time points
        iF  = 1:length(waveletFreqs); % # of freqs.
        iChanSave = 1;

        % chan X event X freq X time
        % make power 10*log(power)
        powerMat(iChanSave,iEv,iF,iT) = 10*log10(rawPow);
        phaseMat(iChanSave,iEv,iF,iT) = rawPhase;

    %     for each eegfile stem, z-score each channel and frequency
        fprintf(' [%.1f sec] --> z-score', toc);  tic;
        stemList = unique({eventTrigger.eegfile});
        
        % indices of the powerMat to Z-score wrt
        fixOnToOff = abs(eventsTriggerXlim(1))*resampledrate:(abs(eventsTriggerXlim(1)) + 1)*resampledrate - 1; % -1 sec to 0 seconds probe word on
        for iStem=1:length(stemList),
            fprintf('.');
            iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));
            for iF = 1:length(waveletFreqs),
%                 allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
                allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,fixOnToOff)),length(iEvStem)*length(fixOnToOff),1); % normalize wrt fixation period
                mu = mean(allVal); stdev = std(allVal);

                % create the power matrix
                powerMatZ(iChanSave,iEvStem,iF,iT) = (powerMat(iChanSave,iEvStem,iF,iT)-mu)/stdev;

                if sum(isnan(powerMatZ(iChanSave,iEvStem,iF,iT)))>0
                    keyboard;
                end
            end
        end
        % set two paramters from robust spectrotemp to 0 
        tWin = 0;
        freq = 0;
        
        fprintf(' [%.1f sec]', toc); tic;
        clear rawPow rawPhase
        disp('powerMatZ, powerMat and phaseMat are created')
    else % OPTION 2: Perform robust spectrotemporal pursuit instead 
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
    clear powerMat
    
    % create vector of the actual seconds in time axis for the powerMat
    % (since its time binned)...
%     LOWERTIME = 1001;
%     UPPERTIME = 6000;
    OVERLAP = 100;
%     FS = 1000;
%     TIMEZERO = 2000;
    if tWin == 0, % if not set yet
        tWin = (LOWERTIME) :OVERLAP/FS: (UPPERTIME);
    end
    timeZero = abs(0-(LOWERTIME))/(OVERLAP/FS);
%     
    % remake powerMatZ to the points that we want (before probe on -> 3.5
    % seconds later
    powerMatZ = squeeze(powerMatZ); % only get the powerMatZ time points we want... (1001 - 2000+3500) = -1.0 seconds -> 3.5 seconds
%     powerMatZ = powerMatZ(:,:,LOWERTIME:UPPERTIME);
    
    if DEBUG,
        size(powerMatZ)
    end
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
    
    %% SPLIT INTO SESSIONS AND BLOCKS
    subjSessions = unique({events.sessionName}); % e.g. sessions 0, 1, 2
    subjBlocks = unique({events.blocknumber});   % e.g. blocks 0,1,2,3,4,5
    
    for iSesh=1:length(subjSessions),
        for iBlock=1:length(subjBlocks),
            sessionBlockIndices = strcmp({events.sessionName}, subjSessions(iSesh)) & ...
                                    strcmp({events.blocknumber}, subjBlocks(iBlock));
                                
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
                    sessionNumber = sessionBlockWordPairEvents(1).sessionNum;
                    
                    thisPowMat = powerMatZ(eventIndices,:,:);
                    
                    %% SAVE PROCESSED DATA IN A MATLAB STRUCT
                    if SAVE,
                        %%- Save this new power matrix Z-scored into data .mat file
                        data.probeWords = THIS_PROBE;                   % the probe words for all events in this struct
                        data.targetWords = THIS_TARGET;                 % the target words for all events in this struct
                        data.sessionNum = sessionNumber;                % the session number
                        data.blockNum = subjBlocks{iBlock};             % the block number
                        data.eegWaveV = eegWaveV(eventIndices,:);                       % eeg wave form
                        data.eegWaveT = eegWaveT;                       % time series for eeg voltage wave
                        data.chanNum = thisChan;                        % store the corresponding channel number
                        data.chanStr = thisChanStr;                     % the string name of the channel
                        data.freqBandYtick = 1:length(freqBandYticks);            % store frequency bands if using wavelet transform
                        data.freqBandYlabel = {freqBandAr.name};
                        data.descriptor = 'Initial processing -2 seconds to 5 seconds after probeWordOn. Time binned with 500ms window and 100ms overlap from -1.1 seconds -> 3.5 seconds';
                        
                        % to plot the axes
%                         set(gca, 'YTick', 1:7, 'YTickLabel', {freqBandAr.name})
                        
                        data.timeZero = timeZero; %ceil((TIMEZERO-LOWERTIME)/OVERLAP);
                        data.vocalization = data.timeZero + ceil([sessionBlockWordPairEvents.responseTime]/OVERLAP);
                        data.powerMatZ = thisPowMat;            % save the condensed power Mat Z-scored
                        data.waveT = tWin;                      % ROBUSTSPECT: save the binned Wave T
                        data.freq = freq;                       % ROBUSTSPECT: save the frequency points

                        %%- SAVING DIR PARAMETERS
                        if ROBUST_SPEC,
                            TYPE_SPECT = 'robust_spec';
                        else
                            TYPE_SPECT = 'morlet_spec';
                        end
                        
                        chanFileName = strcat(num2str(thisChan), '_', thisChanStr, '_', TYPE_SPECT);
                        wordpair_name = strcat(THIS_PROBE, '_', THIS_TARGET);
                        
                        % data directories to save data into
                        homeDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/';
                        homeDir = '/home/adamli/paremap/';
                        dataDir = strcat('condensed_data_', subj);
                        typeTransformDir = fullfile(homeDir, dataDir, TYPE_SPECT);
                        fileDir = fullfile(typeTransformDir, subjSessions{iSesh}, subjBlocks{iBlock}, wordpair_name);
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
    
%     % normalized to 1 so all can be shifted to fit on a single plot
%     eegWaveMeanSub  = eegWaveV-mean(mean(eegWaveV));   %double mean and double max because multiple events from same channel should be normalized together
%     eegWaveShift    = (iChanSave-1)*2 + eegWaveMeanSub./max(max(abs(eegWaveMeanSub)));
%     eegInstPow      = abs(eegWaveMeanSub).^2;
%     eegInstPowShift = (iChanSave-1)*2 + eegInstPow./max(max(abs(eegInstPow)));
%     wavesSft(iChanSave,iEv,iT) = eegWaveShift;
end