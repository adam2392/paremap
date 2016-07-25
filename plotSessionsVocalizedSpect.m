% plot spectrogram per targetword
% put into struct and plot spectrogram

subj = 'NIH034';
typeTransform = 'multitaper';
referenceType = 'bipolar';
blocksComp = 'within_blocks';

addpath('./preprocessing');

subjDataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(subjDataDir, strcat(typeTransform, '_', referenceType, '_', 'vocalization_sessiontargetwords'));
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};
sessions

% load in an example file to get the labels, ticks and timeZero
pairDirs = dir(fullfile(dataDir, sessions{1}, blocks{1}));
exampleDir = fullfile(dataDir, sessions{1}, blocks{1}, pairDirs(4).name);
channelData = dir(exampleDir);
data = load(fullfile(exampleDir, channelData(4).name));
data = data.data;
timeTicks = data.waveT(:,2);

ticks = 1:5:length(timeTicks);
labels = timeTicks(ticks);
timeZero = data.timeZero;
%% array of frequency bands
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
freqBandYtickLabels = {freqBandAr.label};

% eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap';     % work
% eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
% eegRootDirHome = '/Volumes/NIL_PASS';
% eegRootDirJhu = '/home/adamli/paremap';
% 
% % Determine which directory we're working with automatically
% if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
% elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
% elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
% else   error('Neither Work nor Home EEG directories exist! Exiting'); end
% 
% % Either go through all the sessions, or a specific session
% disp(['STEP 1: Going through all sessions with ', typeTransform, ' transform'])
% session = 'Meta Session [all]';
% behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap');
% 
% subjDir = fullfileEEG(eegRootDir,subj); % directory to subject (e.g. NIH034)
% docsDir = fullfileEEG(subjDir,'docs');  % directory to the docs (electordes.m, tagNames.txt, etc.)
% talDir  = fullfileEEG(subjDir,'tal');
% defaultEEGfile = fullfileEEG('/Volumes/Shares/FRNU/data/eeg/',subj,'/eeg.reref/');  % default event eegfile fields point here... switch to local before loading
% 
% %%-Load in the Events For This Task/Patient/Session
% events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
% load(sprintf('%s/events.mat',behDir));  %%- load the events file
% fprintf('Loaded %d events from %s\n', length(events), behDir);
% 
% %%- GET CORRECT EVENTS ONLY
% % POST MODIFY EVENTS based on fields we want (e.g. is it correct or not)?
% correctIndices = find([events.isCorrect]==1);
% events = events(correctIndices);
% [chanList, chanStr, numChannels, eventEEGpath] = loadChannels(docsDir, talDir, THIS_REF_TYPE, USE_CHAN_SUBSET);


%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- SAVING FIGURES OPTIONS
matDir = strcat('./Figures/', subj, '/reinstatement_mat/', typeTransform, '_', ...
    referenceType, '_', blocksComp, '_vocalizationWord/channels_vocalized_session/');

sessionMats = dir(strcat(matDir, '*.mat'));

figure;
for iMat=1:length(sessionMats)
    name = sessionMats(iMat).name
    sessionFile = fullfile(matDir, sessionMats(iMat).name);
    data = load(sessionFile);
    allVocalizedWords = data.allVocalizedWords
    sessionPowerMat = data.sessionPowerMat
    
    %%- loop through number of channels and save in corresponding session
    numChans = size(sessionPowerMat{1}, 2)/7;
    
    seshDir = fullfile(matDir, name(1:end-4));
    if ~exist(seshDir, 'dir')
        mkdir(seshDir);
    end
    
    for iChan=1:numChans
        %%- plot spectrograms and label
        subplot(511);
        powerMat = sessionPowerMat{1};
        featureRange = (iChan-1)*7+1:(iChan-1)*7+7
        imagesc(squeeze(mean(powerMat(:,featureRange,:),1)));
        hold on; colormap(jet); 
        hCbar = colorbar('east');
        title(['Spectrogram for ', ' CLOCK timelocked to vocalization']); %chanStr(iChan),
        xlabel('Time (seconds)');
        set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right')
        % set the heat map settings
        set(gca,'ytick',[1:7],'yticklabel',freqBandYtickLabels)
        set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
        ax = gca;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        
        subplot(512);
        powerMat = sessionPowerMat{2};
        featureRange = (iChan-1)*7+1:(iChan-1)*7+7
        imagesc(squeeze(mean(powerMat(:,featureRange,:),1)));
        hold on; colormap(jet); 
        hCbar = colorbar('east');
        title(['Spectrogram for ', ' JUICE timelocked to vocalization']); %chanStr(iChan),
        xlabel('Time (seconds)');
        set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right')
        % set the heat map settings
        set(gca,'ytick',[1:7],'yticklabel',freqBandYtickLabels)
        set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
        ax = gca;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        
        subplot(513);
        powerMat = sessionPowerMat{3};
        featureRange = (iChan-1)*7+1:(iChan-1)*7+7
        imagesc(squeeze(mean(powerMat(:,featureRange,:),1)));
        hold on; colormap(jet); 
        hCbar = colorbar('east');
        title(['Spectrogram for ', ' PANTS timelocked to vocalization']); %chanStr(iChan),
        xlabel('Time (seconds)');
        set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right')
        % set the heat map settings
        set(gca,'ytick',[1:7],'yticklabel',freqBandYtickLabels)
        set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
        ax = gca;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        
        subplot(514);
        powerMat = sessionPowerMat{4};
        featureRange = (iChan-1)*7+1:(iChan-1)*7+7
        imagesc(squeeze(mean(powerMat(:,featureRange,:),1)));
        hold on; colormap(jet); 
        hCbar = colorbar('east');
        title(['Spectrogram for ', ' BRICK timelocked to vocalization']); %chanStr(iChan),
        xlabel('Time (seconds)');
        set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right')
        % set the heat map settings
        set(gca,'ytick',[1:7],'yticklabel',freqBandYtickLabels)
        set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
        ax = gca;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        
        subplot(515);
        powerMat = sessionPowerMat{5};
        featureRange = (iChan-1)*7+1:(iChan-1)*7+7
        imagesc(squeeze(mean(powerMat(:,featureRange,:),1)));
        hold on; colormap(jet); 
        hCbar = colorbar('east');
        title(['Spectrogram for ', ' GLASS timelocked to vocalization']); %chanStr(iChan),
        xlabel('Time (seconds)');
        set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right')
        % set the heat map settings
        set(gca,'ytick',[1:7],'yticklabel',freqBandYtickLabels)
        set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
        ax = gca;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        
%         %%- SAVE THE FIGURE AFTER CHANGING IT
%         fig = gcf;
%         fig.PaperUnits = 'inches';
%         pos = [0.35, 3.65, 12.55, 7.50];
%         fig.PaperPosition = pos;
%         
%         %%- Save Image
%         figureFile = fullfile(seshDir, chanStr(iChan));
%         print(figureFile, '-dpng', '-r0')
%         savefig(figureFile)
    end
end