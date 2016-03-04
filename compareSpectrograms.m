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

%% Import and load 2 data matrices, or load big matrix with all triggers
if ~exist('data')
    clear all; 
    clc;
    channel_num = '1';
    chanStr = 'G1-global';
    dataDir = 'condensed_data/';
    filename = strcat(dataDir, channel_num, '_', chanStr);

%     data = load('48_MST2-global.mat');
    data = load(filename);
    data = data.data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];

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

correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);

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

figFontAx = 10;

% events x freq x time
powerMatZ = squeeze(powerMatZ);

%% GET TRIGGERS WE WANT
%%- As of 03/04/16: get each trigger word by itself, then word A with B and
%%C and then word B with C
% loop through each trigger and create events struct to pass to plotting
TRIGGER_TYPES   = {'BRICK','CLOCK','GLASS','JUICE','PANTS'}; %-fixation and blockStart not ready for physio analysis
sampEventsMeta = events;  % includes assocaited + and *
probeWords = {sampEventsMeta.probeWord};
targetWords = {sampEventsMeta.targetWord};

for i=1:length(TRIGGER_TYPES)
    THIS_TRIGGER = TRIGGER_TYPES{i};
    
    %%- 01: GET TRIGGER INDICES WE WANT
    switch THIS_TRIGGER,
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
    
    %%- 02: PLOT SPECTROGRAM FOR EACH TRIGGER
    tempFig = figure;
    % loop through each unique trigger for a specific probeword
    for i=1:length(targets);
        % find event indices for this trigger matched with a specific
        % targetword
        targetWord = targets{i};
        eventInd = find(strcmp({sampEventsMeta.probeWord},THIS_TRIGGER) & strcmp({sampEventsMeta.targetWord},targetWord));

        %%- pass in events, powerMatZ, indices of events we want
        %%- Spectrograms of each one
        thisPowMat = powerMatZ(eventInd,:,:);
        powPlot = mean(thisPowMat, 1); % average across events
        titleStr = sprintf('%s EVENTS vs. Target Word(%s) : mean power: chan %s, %d events', THIS_TRIGGER, targetWord, chanStr, size(thisPowMat,1));
        powPlot = squeeze(powPlot);

        % plot spectrogram
        subplot(length(targets), 1, i)
        hImg = imagesc(waveT, 1:7, powPlot);
        hold on; colormap(jet);

        hCbar = colorbar('east');
        % set the heat map settings
        set(gca,'ytick',1:7,'yticklabel',freqBandYtickLabels)
        set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
        set(gca,'fontsize',figFontAx+3)
        set(gca,'XTick',[],'Box','off');
        title(titleStr, 'fontsize',20)
    end
    
    %%- SAVE FIGURE as .mat and as png
    figname = strcat('Figures/', THIS_TRIGGER, '_WithTargets_', chanStr);
    save(figname, 'tempFig');
%     saveas(tempFige, figname, '.fig');
    saveas(tempFig, figname, 'png');
    
    %%- 03: COMPUTE ANOVA?
end

disp('Finished extracting events for the triggers we want!')


%% EXTRACT TRIGGERS
% In this case BRICK, CLOCK, GLASS, PANTS, JUICE
numUniqueTrig = length(unique(trigType));
for thisTrig = 1:numUniqueTrig,
    for otherTrig = thisTrig:numUniqueTrig,
        % only analyze triggers not the same
        if otherTrig ~= thisTrig
            %%- For thisTrig
            % find all event indices for thisTrig (e.g. brickprobe)
            iTrig = find(strcmp(trigType,uniqueTrigType(thisTrig)));
            
            %%- For otherTrig
            % find all event indices for thisTrig (e.g. brickprobe)
            iotherTrig = find(strcmp(trigType,uniqueTrigType(otherTrig)));
            
            % get the two trigger event time series powermatz to test
            firstEvent = newPowerMatZ(iTrig,:,:);
            secondEvent = newPowerMatZ(iotherTrig,:,:);
            
            %% TEST STATISTIC MATRIX
            %%- Loop through newPowerMatZ and generate test statistic per freq/time
            ttestMatZ = zeros(size(newPowerMatZ,2), size(newPowerMatZ,3));
            wstestMatZ = zeros(size(newPowerMatZ,2), size(newPowerMatZ,3));
            for t_index=1:size(newPowerMatZ,3)
                for f_index=1:size(newPowerMatZ,2) 
                    % take all events in this time/freq. slot and compute test statistic
                    [h, p] = ttest2(firstEvent(:,f_index,t_index), secondEvent(:,f_index,t_index));
                    ttestMatZ(f_index,t_index) = p;
                    [p, h] = ranksum(firstEvent(:,f_index,t_index), secondEvent(:,f_index,t_index));
                    wstestMatZ(f_index,t_index) = p;
                end
            end
            
            %%- PLOTTING
            figure
            subplot(2,1,1)
            imagesc(binnedWaveT,log10(waveletFreqs), ttestMatZ)
            hold on;  colormap(jet);
            line([timeZero timeZero], get(gca,'YLim'), 'Color', 'black')
            hCbar = colorbar('east');
            set(hCbar,'ycolor',[1 1 1]*.1, 'fontsize', figFontAx-3, 'YAxisLocation', 'right')
            set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
            set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
            set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
            set(gca,'fontsize',figFontAx+3)
            set(gca,'XTick',[],'Box','off');
            set(gca,'clim',[-3 3]) 
            
            subplot(2,1,2)
            imagesc(binnedWaveT,log10(waveletFreqs),wstestMatZ)
            hold on;  colormap(jet);
            hCbar = colorbar('east');
            set(hCbar,'ycolor',[1 1 1]*.1, 'fontsize', figFontAx-3, 'YAxisLocation', 'right')
            set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
            set(gca,'ytick',log10(freqBandYticks),'yticklabel',freqBandYtickLabels)
            set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
            set(gca,'fontsize',figFontAx+3)
            set(gca,'XTick',[],'Box','off');
            set(gca,'clim',[-3 3]) 
            
            %%- Save files
            firstname = uniqueTrigType{thisTrig};
            secondname = uniqueTrigType{otherTrig};
            if ~exist(strcat(firstname,'_',secondname,'_ttestMat'), 'file') == 2
                save(strcat(firstname,'_',secondname,'_ttestMat'), 'ttestMatZ');
                save(strcat(firstname,'_',secondname,'_wstestMat'), 'wstestMatZ');
            end
        end %end of if
    end 
end
