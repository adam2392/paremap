%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------    VARIABLES USED WHEN RUNNING AS SCRIPT (commment out otherwise)   -----------------------------------%%
% Used to rank cosine similarities based on the features inside feature
% vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% clc
% close all
function rankVocalizationCosineSim(subj, typeTransform, timeLock, referenceType, blocksComp)
% close all
% clc;
% clear all;

addpath('./preprocessing');
% subj = 'NIH034';
% typeTransform = 'multitaper';
% timeLock = 'vocalizationWord';
% referenceType = 'bipolar';
% blocksComp = 'across_blocks';
test = 'postvocal_1';

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

%% load in an example file to get the -> labels, ticks and timeZero
subjDataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(eegRootDir, subjDataDir, strcat(typeTransform, '_', referenceType, '_', 'vocalization_sessiontargetwords'));
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
allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                    'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                    'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', 'GLASS_JUICE', ...
                    'GLASS_GLASS', 'JUICE_JUICE'};

clear chan1 chan2 chanFile chanNums chanRefs eventEEGpath chanTags correctIndices ...
    eventEEGpath eegRootDir eegRootDirHome eegRootDirJhu eegRootDirWork ...
    iChan jackSheet talDir 

%% RUN ACTUAL ANALYSIS
%%- OPEN UP LOG FILE TO PRINT OUT MOST IMPORTANT FEATURES
logFile = strcat(featureMatDir, '/', subj, '.txt');
fid = fopen(logFile, 'w');
LT = 1.5;
newFigDir = fullfile(subjFigDir, 'importantFeatures', strcat(typeTransform, '_', referenceType, '_', timeLock, '_', blocksComp));

allmaxValues = [];

if ~exist(newFigDir)
    mkdir(newFigDir);
end

allVocalizedCosineSim = cell(15,1);

% loop through all session mat files -> extract same, reverse, different
for iMat=1:length(sessionMats),
    sessionMats{iMat}
    data = load(strcat(featureMatDir, '/', sessionMats{iMat}));
    
    %%- each vocalization mat has 15 cells - 1 for each target
    %%word-comparison
    allVocalizedIndices = data.allVocalizedIndices;
%     allVocalizedCosineSim = cell(15,1);
    close all;
    figure;
    
    %% LOOP through target words for this session matrix
    for iTarget=1:length(data.eventReinMat)
        eventReinMat = data.eventReinMat{iTarget};
        featureReinMat = data.featureReinMat{iTarget};
        
        vocalization = zeros(size(featureReinMat, 1), 1);
        for i=1:size(featureReinMat, 1)
            timePeriod = timeZero+3:timeZero+9;
            % get a mini-square time period before the vocalization
            featureCompVocalization = squeeze(featureReinMat(i, timePeriod, timePeriod));

            % compute the mean difference between same pairs and different
            % pairs -> computes a score for this feature and how different it is
            vocalization(i) = mean(featureCompVocalization(:));
            
            %%- test plot to see if the range is correct. CHECK THAT THE
            %%SQUARE ENCLOSES AN AREA WHERE THIS EFFECT IS OCCURRING
            %%POSTVOCALIZATION.
            if i==1,
                subplot(5, 3, iTarget);
                x1 = timePeriod(1);
                x2 = timePeriod(end);
                y1 = timePeriod(1);
                y2 = timePeriod(end);

                imagesc(squeeze(mean(featureReinMat(:, :, :), 1)));
                title(allVocalizedPairs{iTarget})
                colorbar(); colormap('jet');
                axis square; hold on
                xlabel('Time (seconds)'); ylabel('Time (seconds)');
                ax = gca;
                ax.YTick = ticks; ax.YTickLabel = labels;
                ax.XTick = ticks; ax.XTickLabel = labels;
                set(gca,'tickdir','out','YDir','normal'); set(gca, 'box', 'off');
                plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
                plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)

                % draw enclosing square
                x = [x1, x2, x2, x1, x1];
                y = [y1, y1, y2, y2, y1];
                plot(x, y, 'b-', 'LineWidth', 3);
                x = [x1, x2, x2, x1, x1];
                y = [y1, y1, y2, y2, y1];
                plot(x, y, 'b-');  
            end
        end
%         size(allVocalizedCosineSim{iTarget})
%         size(vocalization)
        
%         allVocalizedCosineSim
%         iTarget
%         allVocalizedIndices
        %%- add the vocalization to the master cell
        if allVocalizedIndices(iTarget) == 1 % prevent from catting an empty vector
            if isempty(allVocalizedCosineSim{iTarget}) 
                allVocalizedCosineSim{iTarget} = vocalization;
            else
                try
                    allVocalizedCosineSim{iTarget} = cat(2, allVocalizedCosineSim{iTarget}, vocalization);
                catch
                    disp('error')
                    size(vocalization)
                    size(allVocalizedCosineSim{iTarget})
                    keyboard
                end
            end
        end
    end % loop through targetwords
end % loop through session mats


% save the example reinstatement, so we know where the square was
fig = gcf;
fig.PaperUnits = 'inches';
pos = [0    0.6667   17.5972   10.4028];
fig.PaperPosition = pos;
figureFile = fullfile(newFigDir, strcat(subj, '_exreinstatement', test));

print(figureFile, '-dpng', '-r0')
savefig(figureFile)

%%- loop through all vocalized words, and create a feature importance map
%%and z-scored map
fig = figure(2);
fig2 = gcf;

z = figure(3);
figz = gcf;
for iTarget=1:length(allVocalizedCosineSim)
    subjPostVocalAvge = squeeze(mean(allVocalizedCosineSim{iTarget}, 2));
    Z = zscore(subjPostVocalAvge);
    
    figure(2);
    subplot(5,3, iTarget);
    plot(subjPostVocalAvge);
    hold on;
    title([subj, ' ', allVocalizedPairs{iTarget}]);
    xlabel('Features (channels and 7 frequency bands)');
    ylim([-0.2 1]);
    ylabel('cosine similarity');
    
    figure(3);
    subplot(5, 3, iTarget);
    hist(Z, 20);
    hold on;
    title([subj, ' z-scored ', allVocalizedPairs{iTarget}]);
    xlabel('cosine similarity');
    ylabel('Frequency (n)');
end

%%- Save Image
figure(2);
fig2.PaperUnits = 'inches';
pos = [0    0.6667   17.5972   10.4028];
fig2.PaperPosition = pos;
figureFile = fullfile(newFigDir, strcat(subj, '_featureImportance', test));

print(figureFile, '-dpng', '-r0')
savefig(figureFile)

figure(3);
figz.PaperUnits = 'inches';
pos = [0    0.6667   17.5972   10.4028];
figz.PaperPosition = pos;
figureFile = fullfile(newFigDir, strcat(subj, '_zscored_featureImportance', test));

print(figureFile, '-dpng', '-r0')
savefig(figureFile)
end