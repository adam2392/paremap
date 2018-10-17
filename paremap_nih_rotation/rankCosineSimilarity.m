%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
% Used to rank cosine similarities based on the features inside feature
% vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% clc
% close all
function rankCosineSimilarity(subj, typeTransform, timeLock, referenceType, blocksComp)

subj = 'NIH034';
typeTransform = 'morlet';
timeLock = 'vocalizationWord';
referenceType = 'bipolar';
blocksComp = 'within_blocks';

if strcmp(subj, 'NIH034')
    sessions = {'session_1', 'session_2'};
elseif strcmp(subj, 'NIH039')
    sessions = {'session_0', 'session_1', 'session_3'};
end

eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
eegRootDirHome = '/Volumes/NIL_PASS';
eegRootDirJhu = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

% set the directories with the reinstatements
subjFigDir = fullfile(eegRootDir, 'Figures', subj);
reinstatementDir = strcat(typeTransform, '_', referenceType, '_', blocksComp, '_', timeLock);
featureMatDir = fullfile(subjFigDir, strcat('/reinstatement_mat/', reinstatementDir));

%%- FOR WITHIN BLOCKS RIGHT NOW
%%- GET LIST OF MAT FILES TO EXTRACT DATA FROM
% cd(featureMatDir);
sessionMats = dir(strcat(featureMatDir, '/*.mat'));
sessionMats = {sessionMats.name};

% load in an example file to get the -> labels, ticks and timeZero
THIS_REF_TYPE = referenceType; 
TYPE_TRANSFORM = strcat(typeTransform, '_', referenceType);

% load in an example file to get the -> labels, ticks and timeZero
exDir = fullfile(eegRootDir, 'condensed_data_NIH039/morlet_bipolar/vocalization_allPairs/BRICK_CLOCK/2  3_G2-G3.mat');
[ticks, labels, timeZero] = getPlotMetaData(exDir);

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
USE_CHAN_SUBSET=0;
[chanList, chanStr, numChannels, eventEEGpath] = loadChannels(docsDir, talDir, referenceType, USE_CHAN_SUBSET);

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
newFigDir = fullfile(subjFigDir, 'importantFeatures', strcat(typeTransform, '_', referenceType, '_', timeLock, '_', blocksComp));

allmaxValues = [];

if ~exist(newFigDir)
    mkdir(newFigDir);
end

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
        timePeriod = timeZero-3:timeZero-1;
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

end