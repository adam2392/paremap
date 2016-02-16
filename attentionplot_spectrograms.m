clear all; 
clc;
%%- SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];

% %%- EXTRACTION OPTIONS
% TRIGGER_TYPES   = {'BRICK','CLOCK','GLASS','JUICE','PANTS'}; %-fixation and blockStart not ready for physio analysis
% THIS_TRIGGER    = TRIGGER_TYPES{1};   %%%%%%%%%% SELECT TRIGGER TYPE FOR EXTRACTION AND ANALYSIS
%%- EXTRACTION OPTIONS
TRIGGER_TYPES   = {'sampleStart' 'testStart' 'testEnd' 'asterisk' 'sampleVSastStart' 'sampleVSastEnd' 'fixation' 'blockStart'}; %-fixation and blockStart not ready for physio analysis
THIS_TRIGGER    = TRIGGER_TYPES{1};   %%%%%%%%%% SELECT TRIGGER TYPE FOR EXTRACTION AND ANALYSIS

REF_TYPES       = {'noreref', 'bipolar', 'global'};
THIS_REF_TYPE   = REF_TYPES{3}; % (1) noreref, (2) bipolar, (3) laplacian

HIDE_FIGURES    = 0;
USE_CHAN_SUBSET = 1; %0=all channels (not the subset); >=1 means process than many of the subset
FIG_OFFSET = 0;

%%- FILTERING OPTIONS
BP_FILTER_RAW                 = 1;  %-0 or 1: apply a bandpass filter to the raw traces (1-499 hz)
PROCESS_CHANNELS_SEQUENTIALLY = 1;  %0 or 1:  0 means extract all at once, 1 means sequentially

%%- subset channel I want to extract
iChanListSub  = [48 1];
%%- select all channels, or part of the subset of channels
if USE_CHAN_SUBSET==0,
    iChanList = 1:size(chanList,1);  %all possible channels
else
    iChanList = iChanListSub(1:min([length(iChanListSub) USE_CHAN_SUBSET]));    %select subset of channels
end


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
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/attentionTask/');
    sessStr = '[all]';
else                                  % one session
    disp('STEP 1: Going through one session')
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/attentionTask/', session);
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
%%--------------------------------------------      Create Channel List          -------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        iChanListSub  = [48 1];            %G1, G2, LF1, AST1,
    otherwise
        fprintf('Error, no referencing scheme selected');
end

%%- select all channels, or part of the subset of channels
if USE_CHAN_SUBSET==0,
    iChanList = 1:size(chanList,1);  %all possible channels
else
    iChanList = iChanListSub(1:min([length(iChanListSub) USE_CHAN_SUBSET]));    %select subset of channels
end

%%%%%%%%%%%%%%%%- STEP 2: Manipulate variables for printing and display %%%%%%%%%%%%%%%%%
%%- Get 1. # of channels to use, 
%%-     2. list of channels
%%-     3. list of channel names
chanListUse = [];  chanStrUse = {};
for iChan=iChanList,
    chanListUse(end+1,:) = chanList(iChan,:);
    chanStrUse{end+1}    = chanStr{iChan};
end
% reset variables
chanList = chanListUse; chanStr  = chanStrUse;
clear chanListUse chanStrUse
numChannels = size(chanList,1);
fprintf('\n');
fprintf('STEP 2 -- %d channels to process for %s : %s', numChannels, subj, session);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 3: CREATE EVENT TRIGGERS                  ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampEventsMeta = events;  % includes assocaited + and *
sampEventsTrig = sampEventsMeta(find([sampEventsMeta.isSample])); % get 'issamples'

%*?? the sampeEventsTrig number of elements is different?
% ?? how to find numUniqueTrig? unique Triggers, because that compresses the number of evnets 
% ?? reduce dimensionality of number of events
switch THIS_TRIGGER,
    case 'blockStart'
    case 'sampleStart'
        disp(['Looking at trigger: ', THIS_TRIGGER]);
        eventsMeta = sampEventsMeta;    
        metaYval = [eventsMeta.sampleCount];    
        metaZeroMS = [eventsMeta.msTextStart];  
        metaYstr = 'Sample Onset';     
        metaYtransition = find(diff([eventsMeta.asteriskType]));
        eventsTrig = sampEventsTrig;    
        trigType = [eventsTrig.asteriskType];   
        trigZeroMS = [eventsTrig.msTextStart];  trigTypeStr = 'asteriskType';
        eventsTriggerXlim = [-2.25 5.25];
        eventsAveWindowMS = [-1000 -500; -500 0; 0 500; 500 1000; 1000 1500; 1500 2000; 2000 2500; 3500 4000]; % list of time windows over which EEG data is averaged for t-tests
        %eventsTriggerXlim = [min(min(eventsAveWindowMS)) max(max(eventsAveWindowMS))]/1000;
    case 'testStart'
    case 'testEnd'
    case 'fixation'
    case 'asterisk'
    case 'sampleVSastStart'
    case 'sampleVSastEnd'
        
    otherwise
        error('no event trigger selected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 4: DATA INPUT TO GETE_MS, MULTIPHASEVEC3 AND ZSCORE   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%- Input to gete_ms
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

% remap event pointer from default (server) to local copy of the EEG data
for iEvent=1:length(eventTrigger),
    eventTrigger(iEvent).eegfile = regexprep(eventTrigger(iEvent).eegfile,defaultEEGfile,fullfileEEG(subjDir,eventEEGpath));
end

%%- gets the range of frequencies using eeganalparams
waveletFreqs = eeganalparams('freqs');
waveletWidth = eeganalparams('width');

%-- pre-allocate memory (to make sure it can be done!)
if PROCESS_CHANNELS_SEQUENTIALLY==0,  numChanPrealloc = numChannels;  else  numChanPrealloc = 1; end
% #channels X #events X #freqs. X #timepoints = 4D array
arrayGB = numChanPrealloc * length(eventTrigger) * length(waveletFreqs) * DurationMS * 8 / 2^30 % 8 bytes per double, 2^30 bytes/GB
% initialize power matrices to make sure they can be stored
powerMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
powerMatZ = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
phaseMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 5: Loop through the channels: extract, filter, processes, and save...   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numUniqueTrig = 0; % to just test for 2 plots for now

% loop through channels
for iChan=1:numChannels 
    thisChan = chanList(iChan,:);   % the channel to use in this loop (e.g. 48)
    
    %%- gete_ms: get the eegWaveV
    eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);

    % notch filter to eliminate 60 Hz noise
    fprintf(' [%.1f sec] --> notch filt', toc); tic;
    eegWaveV = buttfilt(eegWaveV,[59.5 60.5],resampledrate,'stop',1); %-filter is overkill: order 1 --> 25 dB drop (removing 5-15dB peak)
    
    %%- multiphasevec3: get the phase and power
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
    
    size(eegWaveT)
    size(eegWaveV)
    size(rawPow)
    size(rawPhase)
    
    % temp indicies
    iEv = 1:length(eventTrigger); % # of events
    iT  = 1:size(eegWaveV,2); % # of time points
    iF  = 1:length(waveletFreqs); % # of freqs.
    iChanSave = 1;
    
    % chan X event X freq X time
    powerMat(iChanSave,iEv,iF,iT) = 10*log10(rawPow);
    phaseMat(iChanSave,iEv,iF,iT) = rawPhase;
    
    size(powerMat)
    size(phaseMat)
    
    % for each eegfile stem, z-score each channel and frequency
    fprintf(' [%.1f sec] --> z-score', toc);  tic;
    stemList = unique({eventTrigger.eegfile});
    for iStem=1:length(stemList),
        fprintf('.');
        iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));
        for iF = 1:length(waveletFreqs),
            allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
            mu = mean(allVal); stdev = std(allVal);
            powerMatZ(iChanSave,iEvStem,iF,iT) = (powerMat(iChanSave,iEvStem,iF,iT)-mu)/stdev;
            if sum(isnan(powerMatZ(iChanSave,iEvStem,iF,iT)))>0
                keyboard;
            end
        end
    end
    fprintf(' [%.1f sec]', toc); tic;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%------------------ STEP 5i: Plot Evoked   ---------------------------------------%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iTrig = size(powerMatZ,2); % just get all event triggers instead of unique ones
    
    % normalized to 1 so all can be shifted to fit on a single plot
    eegWaveMeanSub  = eegWaveV-mean(mean(eegWaveV));   %double mean and double max because multiple events from same channel should be normalized together
    eegWaveShift    = (iChanSave-1)*2 + eegWaveMeanSub./max(max(abs(eegWaveMeanSub)));
    eegInstPow      = abs(eegWaveMeanSub).^2;
    eegInstPowShift = (iChanSave-1)*2 + eegInstPow./max(max(abs(eegInstPow)));
    
    wavesSft(iChanSave,iEv,iT) = eegWaveShift;
    
    % x-axis of time series
    waveT = eegWaveT;
    
    %-compute the evoked response for each trigger type
    wavePlot = wavesSft(:,iTrig,:);
    clear evokedRespMu evokedRespSEM;
    chanNum = iChan;
    evokedRespMu(1:size(wavePlot,3))  = mean(wavePlot(chanNum,:,:),2);
    evokedRespSEM(1:size(wavePlot,3)) = std(wavePlot(chanNum,:,:),0,2)./sqrt(size(wavePlot,2));
    
    cTrigAr = 'rkbgmcy';
    thisTrig = 1;
    cTrig = cTrigAr(thisTrig);
    
    % actually plot evoked mean
    subplot(1+numUniqueTrig,1,1)
    axWave = gca;
    hEvoked = plot(waveT, evokedRespMu, 'k-'); hold on
    set(hEvoked,'Color',cTrig,'LineWidth',2)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%------------------ STEP 5ii: Plot Spectrograms   ---------------------------------------%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thisPowMat = powerMatZ(chanNum,1:iTrig,:,:); 
    if thisTrig==1, 
        cAx = [-.8 .8]; 
        if length(iTrig)<4,
            cAx=cAx*3; 
        elseif length(iTrig)>100, 
            cAx=cAx/2; 
        end;   
    end;   % z-scored power... set cAx for first panel and keep same for all others
    
    %-array of frequency bands
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

    for iFB=1:length(freqBandAr),
        freqBandAr(iFB).centerF = mean(freqBandAr(iFB).rangeF);
        %freqBandAr(iFB).label   = sprintf('%s-%.0fHz', freqBandAr(iFB).name(1:[ min( [length(freqBandAr(iFB).name), 6] )]), freqBandAr(iFB).centerF);
        freqBandAr(iFB).label   = sprintf('%s [%.0f-%.0f Hz]', freqBandAr(iFB).name, freqBandAr(iFB).rangeF);
    end

    freqBandYticks  = unique([freqBandAr(1:7).rangeF]);
    for iFB=1:length(freqBandYticks), freqBandYtickLabels{iFB} = sprintf('%.0f Hz', freqBandYticks(iFB)); end
    
    % reduce the dimensionality of the matrix for plotting
    powPlot = mean(thisPowMat,2); 
    titleStr = sprintf('mean power: chan %s, %d events', chanStr{chanNum}, size(thisPowMat,2));
    powPlot = squeeze(powPlot); % squeeze out the singleton dimension
    
    % actually plot the spectrogram
    subplot(1+numUniqueTrig,1,2)
    axSpec(thisTrig)=gca;
    hImg    = imagesc(waveT,log10(waveletFreqs),powPlot); 
    hold on;  colormap(jet);
    
    % set the heat map settings
    set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
    set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
    set(gca,'fontsize',14)
    set(gca,'XTick',[],'Box','off');
    set(gca,'clim',[-1 1]*0.5)
end