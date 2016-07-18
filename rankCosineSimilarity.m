%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
% Used to rank cosine similarities based on the features inside feature
% vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all

subj = 'NIH034';
typeTransform = 'morlet';
timeLock = 'vocalization';
referenceType = 'bipolar';
typeReinstatement = 'within_blocks';

% set the directories with the reinstatements
subjFigDir = fullfile('Figures', subj);
reinstatementDir = strcat(typeTransform, '_', referenceType, '/', typeReinstatement, '_', timeLock);
featureMatDir = strcat('./Figures/', subj, '/reinstatement_mat/', reinstatementDir);

%%- FOR WITHIN BLOCKS RIGHT NOW
%%- GET LIST OF MAT FILES TO EXTRACT DATA FROM
% cd(featureMatDir);
sessionMats = dir(strcat(featureMatDir, '/*.mat'));
sessionMats = {sessionMats.name};

% load in an example file to get the -> labels, ticks and timeZero
THIS_REF_TYPE = referenceType; 
TYPE_TRANSFORM = strcat(typeTransform, '_', referenceType);
CUE_LOCK = strcat(timeLock);

dataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(dataDir, TYPE_TRANSFORM, CUE_LOCK)
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};
pairDirs = dir(fullfile(dataDir, sessions{1}, blocks{1}));
exampleDir = fullfile(dataDir, sessions{1}, blocks{1}, pairDirs(4).name);
channelData = dir(exampleDir);
data = load(fullfile(exampleDir, channelData(4).name));
data = data.data;
timeTicks = data.waveT(:,2);

ticks = 1:5:length(timeTicks);
labels = timeTicks(ticks);
timeZero = data.timeZero;

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

jackSheet = fullfileEEG(docsDir, 'jacksheetMaster.txt');
[chanNums chanTags] = textread(jackSheet,'%d%s%*s');

%%% always look at all electrodes... worry about "good" and "bad" later (bad means inter-ictal activity or seizure activity)
%- three referencing options:  noreref (should manually subtract reference channel), reref bioploar, and reref laplacian
chanStr = {};   % cell for all the channel names
chanFile = 0;   % file for the channels (e.g. ~/NIH034/tal/leads.txt) 
chanList = [];  % list of the channels (e.g. 1-96)

switch referenceType
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

%%- get the information about each entry in the feature's row (1st
%%dimension), which contains frequency bands and electrodes
freqBands = {'delta', 'theta', 'alpha', 'beta', 'low gamma', 'high gamma', 'HFO'};

clear chan1 chan2 chanFile chanNums chanRefs eventEEGpath chanTags correctIndices ...
    eventEEGpath eegRootDir eegRootDirHome eegRootDirJhu eegRootDirWork ...
    iChan jackSheet talDir 

%%- OPEN UP LOG FILE TO PRINT OUT MOST IMPORTANT FEATURES
logFile = strcat(featureMatDir, '/', subj, '.txt');
fid = fopen(logFile, 'w');
LT = 1.5;
newFigDir = fullfile(subjFigDir, 'importantFeatures', strcat(typeTransform, referenceType, '_', typeReinstatement));

allmaxValues = [];

if ~exist(newFigDir)
    mkdir(newFigDir);
end

sessions = {'session_1', 'session_2'};

% loop through all session mat files -> extract same, reverse, different
for iMat=1:length(sessionMats),
    sessionMats{iMat}
    data = load(strcat(featureMatDir, '/', sessionMats{iMat}));
    featureSame = data.featureSame;
    featureDiff = data.featureDiff;
    
    %%- 01: PRODUCE RANKING METRIC FOR EACH SESSION-BLOCK REINSTATEMENT
    % get the time period before timeZero for this metric and evaluate the
    % average
    preVocal = zeros(size(featureSame, 1), 1);
    for i=1:size(featureSame, 1)
        timePeriod = timeZero-2:timeZero-1;
        % get a mini-square time period before the vocalization
        featureSamePreVocal = squeeze(featureSame(i, timePeriod, timePeriod));
        featureDiffPreVocal = squeeze(featureDiff(i, timePeriod, timePeriod));
        
%         size(featureSamePreVocal)
        
        % compute the mean difference between same pairs and different
        % pairs -> computes a score for this feature and how different it is
        featurePreVocal = featureSamePreVocal - featureDiffPreVocal;
        preVocal(i) = mean(featurePreVocal(:));
    end
    
    % evaluate a single point - w/ multitaper would be good
    
    %%- 02: USING RANKING METRIC, PRODUCT OUTPUT .TXT FILE ON TOP 10
    %%FEATURES CONTRIBUTING TO DIFFERENCES IN SAME/DIFFERENT
    %%PRE-VOCALIZATION
    size(preVocal)
    
    N = 30; % the first N largest feature differences
    % find max indice and the related channel/freq band
    [sortedX, sortingIndices] = sort(preVocal, 'descend'); % sorted in descending order
    maxValues = sortedX(1:N);
    maxValueIndices = sortingIndices(1:N);
    
    if isempty(allmaxValues)
        allmaxValues = preVocal;
    else
        allmaxValues = cat(2, allmaxValues, preVocal);
    end
    %%- plot
%     figure
%     plot(preVocal)
%     title(['Feature Ranking of Cosine Similarity In a PreVocalization Block'])
%     hold on
%     xlabel('features (7 frequency bands for every electrode)');
%     ylabel('Different minus Same Pair Comparison');
%     ax = gca;
    
    %%- plot the preVocal max
%     figure
%     plot(maxValues)

    importantChannelIndices = ceil(maxValueIndices/7);
    importantChannels = chanStr(importantChannelIndices);
    importantFreqIndices = mod(maxValueIndices, 7);
    importantFreqIndices = importantFreqIndices + 1;
    importantFreqs = freqBands(importantFreqIndices);

    fprintf(fid, '\n %6s \n', sessionMats{iMat}); % print which session it came from
    for i=1:N
        fprintf(fid, '%6s \n', [importantChannels{i}, ' , ', importantFreqs{i}]);
    end   
    
    %%- save the important feature indices in a mat file
    toSave{iMat} = [importantChannelIndices, importantFreqIndices];
    
    %%- log any overlapping indices
%     figure
%     fig = {};
%     fig{end+1} = subplot(311)
%     imagesc(squeeze(mean(featureSame(maxValueIndices, :, :),1)));
%     title(['Same Pairs Cosine Similarity for ', sessionMats{iMat}])
%     colorbar();
%     clim = get(gca, 'clim');
%     axis square
%     hold on
%     xlabel('Time (seconds)');
%     ylabel('Time (seconds)');
%     ax = gca;
%     ax.YTick = ticks;
%     ax.YTickLabel = labels;
%     ax.XTick = ticks;
%     ax.XTickLabel = labels;
%     colormap('jet');
%     set(gca,'tickdir','out','YDir','normal');
%     set(gca, 'box', 'off');
%     plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%     plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
%     
%     fig{end+1} = subplot(312);
%     imagesc(squeeze(mean(featureDiff(maxValueIndices, :, :),1)));
%     title(['Different Word Pairs Cosine Similarity for ', sessionMats{iMat}])
%     colorbar();
%     set(gca, 'clim', clim);
%     axis square
%     hold on
%     xlabel('Time (seconds)');
%     ylabel('Time (seconds)');
%     ax = gca;
%     ax.YTick = ticks;
%     ax.YTickLabel = labels;
%     ax.XTick = ticks;
%     ax.XTickLabel = labels;
%     colormap('jet');
%     set(gca,'tickdir','out','YDir','normal');
%     set(gca, 'box', 'off');
%     plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%     plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
% 
%     fig{end+1} = subplot(313);
%     imagesc(squeeze(mean(featureSame(maxValueIndices, :, :),1)) - squeeze(mean(featureDiff(maxValueIndices, :, :),1)));
%     title({'Within-Blocks', ['Same-Different Word Pairs Cosine Similarity']})
%     colorbar();
%     axis square
%     hold on
%     xlabel('Time (seconds)');
%     ylabel('Time (seconds)');
%     ax = gca;
%     ax.YTick = ticks;
%     ax.YTickLabel = labels;
%     ax.XTick = ticks;
%     ax.XTickLabel = labels;
%     colormap('jet');
%     set(gca,'tickdir','out','YDir','normal');
%     set(gca, 'box', 'off');
%     plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%     plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
%     
%     %%- Save Image
%     figureFile = fullfile(newFigDir, sessionMats{iMat}(1:end-4));
%     print(figureFile, '-dpng', '-r0')
%     savefig(figureFile)
end

subjPreVocalAvge = squeeze(mean(allmaxValues, 2));
Z = zscore(subjPreVocalAvge);

size(subjPreVocalAvge)
figure;
subplot(211);
plot(subjPreVocalAvge);
hold on;
title([subj, ' Same-Different Pair Reinstatement Block PreVocalization']);
xlabel('Features (channels and 7 frequency bands)');
ylabel('Reinstatement Difference (metric = cosine similarity)');
subplot(212);
hist(Z, 20);
title([subj, ' ZScored Same-Different Pair Reinstatement Block PreVocalization']);
xlabel('Reinstatement Difference (metric = cosine similarity)');
ylabel('Frequency (n)');

%%- Save Image
fig = gcf;
fig.PaperUnits = 'inches';
pos = [0    0.6667   17.5972   10.4028];
fig.PaperPosition = pos;
figureFile = fullfile(newFigDir, strcat(subj, '_featureImportance'));
print(figureFile, '-dpng', '-r0')
savefig(figureFile)

% save mat file
saveIndicesFile = fullfile(newFigDir, strcat(subj, '_importantIndices'));
save(saveIndicesFile, 'toSave', 'Z', 'allmaxValues');

fclose(fid);