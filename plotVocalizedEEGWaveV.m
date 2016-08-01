subj = 'NIH034';
typeTransform = 'morlet';
referenceType = 'bipolar';
blocksComp = 'within_blocks';

addpath('./preprocessing');
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

eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
eegRootDirHome = '/Volumes/NIL_PASS';
eegRootDirJhu = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

% Either go through all the sessions, or a specific session
disp(['STEP 1: Going through all sessions with ', typeTransform, ' transform'])
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
USE_CHAN_SUBSET=0;
[chanList, chanStr, numChannels, eventEEGpath] = loadChannels(docsDir, talDir, referenceType, USE_CHAN_SUBSET);

subjDataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(eegRootDir, subjDataDir, strcat(typeTransform, '_', referenceType, '_targetWords'));
targetWords = dir(dataDir);
targetWords = {targetWords(3:end).name};
%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%- SAVING FIGURES OPTIONS®
% matDir = fullfile(eegRootDir, strcat('./Figures/', subj, '/spectrograms/', typeTransform, '_', ...
%     referenceType, '_', blocksComp, '_vocalizationWord/channels_vocalized_session/'));
matDir = fullfile(eegRootDir, strcat('./Figures/', subj, '/spectrograms/', typeTransform, '_', ...
    referenceType, '_targetWords/'));
eegWaveDir = fullfile(eegRootDir, strcat('./Figures/', subj, '/eegwaves/', typeTransform, '_', ...
    referenceType, '_targetWords/'));
if ~exist(eegWaveDir, 'dir')
    mkdir(eegWaveDir);
end
numChannels = length(dir(fullfile(dataDir, targetWords{1}, '*.mat')));

LT = 1.5;

for iChan=1:numChannels
    figure;
    FA={};
    clim = [3 -4];
    for iTarget=1:length(targetWords)
        chanFiles = dir(fullfile(dataDir, targetWords{1}, '*.mat'));
        chanFiles = {chanFiles.name};
        chanFile = fullfile(dataDir, targetWords{iTarget}, chanFiles{iChan});
        data = load(chanFile);
        data = data.data;
        targetEEG = data.eegWaveV;
        timeTicks = data.waveT(:,2);

        ticks = 1:1000:length(targetEEG);
        labels = [-2, -1, 0, 1, 2, 3];
        timeZero = 2000;
        
        % plot
        FA{iTarget} = subplot(3,2, iTarget);
        randindice = randsample(size(targetEEG,1), 1);
        plot(squeeze(targetEEG(randindice,:)));
        hold on;
        title(['EEG Voltage Trace for ', chanStr(iChan), ' ', targetWords{iTarget}, ' timelocked to vocalization']); %chanStr(iChan),
        xlabel('Time (seconds)');
        % spectrogram should have low freq on the bottom
        ax = gca;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
    end

    %%- SAVE THE FIGURE AFTER CHANGING IT
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.Units = 'inches';
    pos = [11.9375   -7.4896   19.3021   11.3750];
    fig.Position = pos;
    fig.PaperPosition = pos;

    %%- Save Image
    figureFile = fullfile(eegWaveDir, strcat(chanStr(iChan), '.png'));
    print(figureFile{:}, '-dpng', '-r0')

    pause(0.05);
    close all;
end