%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSUMING DATA IS PASSED IN %%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% PARAMETERS FOR RUNNING PREPROCESS
subj = 'NIH034';
sessNum = [0, 1, 2];

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
PROCESS_CHANNELS_SEQUENTIALLY = 0;  %0 or 1:  0 means extract all at once, 1 means sequentially

%% LOAD EVENTS STRUCT AND SET DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
eegRootDirHome = '/home/adamli/paremap';

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
eventsTriggerXlim = [-1 5]; % range of time to get data from (-2 seconds to 5 seconds after mstime (probeWordOn)) 
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
arrayGB = numChanPrealloc * length(eventTrigger) * DurationMS * 8 / 2^30; % 8 bytes per double, 2^30 bytes/GB

clear defaultEEGfile subjDir eventEEGpath
% print statements for debugging and process checking
fprintf('\n');
fprintf('STEP 4 -- %d events to process for %s : %s', length(eventTrigger), subj, session);
fprintf('\n');
fprintf('The amount of RAM (GB) needed is: %d', arrayGB);
fprintf('\n\n');

fprintf('Number of preallocated channels are: %d', numChanPrealloc)
fprintf('\n');
fprintf('Duration of analysis is: %d\n', DurationMS)

%% LOOP THROUGH EVERY EVENT AND GET A SLIDING CONNECTIVITY MATRIX
tic;
for iEvent=1:length(eventTrigger),
    event = eventTrigger(iEvent);
    
    eegWave = zeros(numChanPrealloc, DurationMS);
    %%- LOOP THROUGH ALL CHANNELS FOR SUBJECT TO GET DATA
    for thisChan=1:numChannels,
        %%- gete_ms: get the eegWaveV
        % eegwaveform for each event over the duration of time for a certain channel
        eegWaveV = gete_ms(thisChan,event,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);

        %%- 01: Apply notch filter to data
        % notch/stop filter at 60 Hz 
        eegWaveV = buttfilt(eegWaveV, [59.5, 60.5], resampledrate, 'stop', 1); 
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area

        %%- 02: Z-score all the data per event (subtract mean and divide by standard
        %%deviation)
        avgeEEG =  mean(eegWaveV, 2);
        stdEEG = std(eegWaveV, 1, 2);
        eegWaveV = bsxfun(@minus, eegWaveV, avgeEEG)./stdEEG(iEvent);
        eegWaveV(isnan(eegWaveV)) = 0;
        eegWave(thisChan, :) = eegWaveV; % store data into channelXtime matrix
    end
    fprintf('Creating channels EEG matrix: %.3s', toc);

    % set the # of points of the fft
    Fs = 1000;
    nfft = 500;
    % set the step size for the sub-sections
    stepsize = round(nfft/5); % overlap 25%
    
    %%- W(j) in Welch's method
    % set window function for periodograms from Welch's method
    convwin = 0.54-0.46 * cos((2*pi*(1:nfft))/(nfft-1));
    convwin = convwin(:); %vectorize this 
    
    
    % set range of frequencies for periodigrams
    frq = (Fs/2) .* linspace(0, 1, nfft/2);
    
    %%- 03: split the data into overlapping sub-sections for cross-power
    % zero-lag cross-correlation for each pair of iEEG electrodes
    nsections = ceil((size(eegWave, 2) - (nfft-stepsize))/stepsize);
    
    %%- 04: Compute cross-power spectra
    CrossPwr = zeros(nfft/2, numChannels*(numChannels+1)/2);
    CrossSpc = zeros(nfft/2, numChannels*(numChannels+1)/2);
    for i=1:nsections,
        % ... extract the current section of data
        datasec = eegWave(:, (i-1)*stepsize+1 : min((i-1)*stepsize+nfft, size(eegWave,2)))';
        
        % ... check if data is all zeros
        if (sum(sum(abs(datasec)))>0),
            % ... comput the fft. size(X) = (nfft/2) x # of channels
            X = fft(repmat(convwin(1 : min(nfft, size(datasec, 1))), 1, numChannels) .* datasec, nfft);
            X = X(1:nfft/2,:);
            size(X)
            
            %%%%%% Read up on Welch's method and then implement...
        end
        
        stepsize
        nsections
        size(datasec)
        
    end
    
end
%     fileDir = strcat('/home/adamli/paremap/', subj, '_eegwaves');
%     subjeegfile = fullfile(fileDir, 'eegData');
%     if ~exist(fileDir)
%         mkdir(fileDir);
%     end
%     save(subjeegfile, 'eegWave', '-v7.3');      

   























