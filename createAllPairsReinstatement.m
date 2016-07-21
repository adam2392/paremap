% Script to run reinstatement construction for pairs of words that are
% vocalized. 
% 1) Run without regard for session separation
% 2) Only compare everything within a session and block

close all;
clc;
clear all;

addpath('./reinstatement_allPairs/');
%% 1) Without session separation
subj = 'NIH034';
timeLock = 'vocalization_allPairs';
referenceType = 'bipolar';
typeTransform = 'morlet';

THIS_REF_TYPE = referenceType; 
TYPE_TRANSFORM = strcat(typeTransform, '_', referenceType);
CUE_LOCK = strcat(timeLock);

[allWordPairs, dataDir] = loadAllPairs(subj, TYPE_TRANSFORM, CUE_LOCK);
wordDirs = fullfile(dataDir, allWordPairs);

%%- SAVING FIGURES OPTIONS
if strcmp(CUE_LOCK, 'vocalization_allPairs')
    figureDir = strcat('./Figures/', subj, '/reinstatement/', TYPE_TRANSFORM,'/within_blocks_vocalization_allPairs/');
    matDir = strcat('./Figures/', subj, '/reinstatement_mat/', TYPE_TRANSFORM,'/within_blocks_vocalization_allPairs/');
end

% make directory, if doens't already exist
if ~exist(figureDir) mkdir(figureDir); end
if ~exist(matDir)    mkdir(matDir);    end
% set linethickness
LT = 1.5;

% load in an example file to get the -> labels, ticks and timeZero
exDir = strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data_', subj, '/morlet_bipolar/vocalization_allPairs/BRICK_CLOCK/2  3_G2-G3.mat');
[ticks, labels, timeZero] = getPlotMetaData(exDir);

%%- 1. LOOP THROUGH EACH WORDPAIR DIRECTORY
for iPair=1:length(allWordPairs)
    wordPairDir = wordDirs{iPair}; % directory to current word Pair
    % get the channel names
    channels = dir(wordPairDir);
    channels = {channels(3:end).name};
    
    %%- A. CREATE SAME PAIR REINSTATEMENT INPUT
    [samePairFeatureMat1, samePairFeatureMat2] = createSamePairReinInput(channels, wordPairDir);
    
    %%- B. CREATE DIFF PAIR REINSTATEMENT INPUT
    wordPair = allWordPairs{iPair};
    allowedIndices = [];
    words = strsplit(wordPair, '_');
    for i=1:length(allWordPairs)
        if (isempty(strfind(allWordPairs{i}, words{1})) && isempty(strfind(allWordPairs{i}, words{2})))
            
            allowedIndices = [allowedIndices; i];
        end
    end
    allWordPairs{allowedIndices};
    randWordPair = randsample(allowedIndices, 1);

    wordPairDir1 = wordPairDir;
    wordPairDir2 = wordDirs{randWordPair};
    [diffPairFeatureMat1, diffPairFeatureMat2] = createDiffPairReinInput(channels, wordPairDir1, wordPairDir2);
    
    %%- C. COMPUTE REINSTATEMENT MAPS
    % get random indices to downsize the comparisons made between pair events
    randIndices = randsample(size(samePairFeatureMat1, 1), min(size(diffPairFeatureMat1, 1), size(samePairFeatureMat2, 1)));
    samePairFeatureMat1 = samePairFeatureMat1(randIndices, :, :);
    samePairFeatureMat2 = samePairFeatureMat2(randIndices, :, :);
    
    randIndices = randsample(size(diffPairFeatureMat1, 1), min(size(diffPairFeatureMat1, 1), size(samePairFeatureMat2, 1)));
    diffPairFeatureMat1 = diffPairFeatureMat1(randIndices, :, :);
    diffPairFeatureMat2 = diffPairFeatureMat2(randIndices, :, :);
    
    size(samePairFeatureMat1)
    size(diffPairFeatureMat1)
    
    fprintf('Computing reinstatement matrices \n');
    [eventSame, featureSame] = compute_reinstatement(samePairFeatureMat1, samePairFeatureMat2);
    [eventDiff, featureDiff] = compute_reinstatement(diffPairFeatureMat1, diffPairFeatureMat2);
    
    %%- D. PLOT REINSTATEMENTS
    figure;
    subplot(311);
    imagesc(squeeze(mean(eventSame(:, :, :),1)));
    title([wordPair, ' Cosine Similarity ', ...
        ' with ', num2str(size(eventSame,1)), ' events'])
    colorbar(); colormap('jet');
    clim = get(gca, 'clim');
    axis square; hold on;
    xlabel('Time (seconds)');
    ylabel('Time (seconds)');
    ax = gca;
    ax.YTick = ticks; ax.YTickLabel = labels;
    ax.XTick = ticks; ax.XTickLabel = labels;
    set(gca,'tickdir','out','YDir','normal', 'box', 'off');
    plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
    plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
    
    subplot(312);
    imagesc(squeeze(mean(eventDiff(:, :, :),1)));
    title([wordPair ' vs. ', allWordPairs{randWordPair}, ' Cosine Similarity ', ...
        ' with ', num2str(size(eventDiff,1)), ' events'])
    colorbar(); colormap('jet');
    clim = get(gca, 'clim');
    axis square; hold on;
    xlabel('Time (seconds)');
    ylabel('Time (seconds)');
    ax = gca;
    ax.YTick = ticks; ax.YTickLabel = labels;
    ax.XTick = ticks; ax.XTickLabel = labels;
    set(gca,'tickdir','out','YDir','normal', 'box', 'off');
    plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
    plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
    
    subplot(313);
    imagesc(squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(eventDiff(:, :, :),1)));
    title(['Same - Different Cosine Similarity ', subj])
    colorbar(); colormap('jet');
    clim = get(gca, 'clim');
    axis square; hold on;
    xlabel('Time (seconds)');
    ylabel('Time (seconds)');
    ax = gca;
    ax.YTick = ticks; ax.YTickLabel = labels;
    ax.XTick = ticks; ax.XTickLabel = labels;
    set(gca,'tickdir','out','YDir','normal', 'box', 'off');
    plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
    plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
    
    %%- Save Image and save reinstatement matrices
    figureFile = strcat(figureDir, allWordPairs{iPair}, 'vs', allWordPairs{randWordPair});
    matFile = strcat(matDir, allWordPairs{iPair});
    comparedWord = allWordPairs{randWordPair};
    save(strcat(matFile, '.mat'), 'eventSame', 'featureSame', ...
                                        'eventDiff', 'featureDiff',...
                                    'comparedWord');
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    pos = [0.35, 3.65, 12.55, 7.50];
    fig.PaperPosition = pos;
    print(figureFile, '-dpng', '-r0')
    savefig(figureFile)

    pause(0.1); 
end

