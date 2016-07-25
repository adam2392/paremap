%         -- Creating matrices that are eventsXtimeXfeatures, which can be
%         fed into compute_reinstatement.m
%         -- This is done for ACROSS blocks analysis of the paremap task

function createAcrossBlocksReinstatementMat(subj, typeTransform, timeLock, referenceType)
close all;
clc;

%% PARAMETERS FOR RUNNING PREPROCESS
% subj = 'NIH034';
% timeLock = 'vocalization';
% referenceType = 'bipolar';
% % winSize = 500;
% % stepSize = 100;
% typeTransform = 'morlet';

expected_timeLocks = {'vocalization', 'matchword', 'probeword'};
expected_transforms = {'morlet', 'multitaper'};
REF_TYPES = {'noreref', 'bipolar', 'global'};
if ~ismember(timeLock, expected_timeLocks)
    disp('timeLock should be vocalization, matchword, or probeword');
end
if ~ismember(referenceType, REF_TYPES)
    disp('reference types are noreref, bipolar, or global');
end
if ~ismember(typeTransform, expected_transforms)
    disp('transform types are morlet, multitaper');
end
THIS_REF_TYPE = referenceType; 
TYPE_TRANSFORM = strcat(typeTransform, '_', referenceType);
CUE_LOCK = strcat(timeLock);

addpath('./m_reinstatement/');

%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjDir = strcat('./condensed_data_', subj);
dataDir = fullfile(subjDir, strcat(typeTransform, '_', referenceType, '_', timeLock));
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

sessions % display which sessions we're working with

%%- SAVING FIGURES OPTIONS
if strcmp(CUE_LOCK, 'vocalization')
    figureDir = strcat('./Figures/', subj, '/reinstatement/', typeTransform, '_', referenceType, '_', 'across_blocks_vocalization/');
    matDir = strcat('./Figures/', subj, '/reinstatement_mat/', typeTransform, '_', referenceType, '_', 'across_blocks_vocalization/');
elseif strcmp(CUE_LOCK, 'matchword')
    figureDir = strcat('./Figures/', subj, '/reinstatement/', typeTransform, '_', referenceType, '_', 'across_blocks_matchword/');
    matDir =strcat('./Figures/', subj, '/reinstatement_mat/', typeTransform, '_', referenceType, '_', 'across_blocks_matchword/');
elseif strcmp(CUE_LOCK, 'probeword')
    figureDir = strcat('./Figures/', subj, '/reinstatement/', typeTransform, '_', referenceType, '_', 'across_blocks_probeon/');
    matDir = strcat('./Figures/', subj, '/reinstatement_mat/', typeTransform, '_', referenceType, '_', 'across_blocks_probeon/');
end

%%- SAVING FIGURES OPTIONS
if ~exist(figureDir, 'dir') mkdir(figureDir); end
if ~exist(matDir, 'dir')    mkdir(matDir);    end

% set linethickness
LT = 1.5;
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

%%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:length(blocks)-1,
        fprintf('%6s \n', strcat('On session ', num2str(iSesh), ' and block ', num2str(iBlock)));
        
        % get word pairs in this session-block (i)
        firstwordpairs = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        firstwordpairs = {firstwordpairs(3:end).name};
        % get the word pairs in the next block (i+1)
        secondwordpairs = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock+1}));
        secondwordpairs = {secondwordpairs(3:end).name};
        
        %%- CREATE ACROSS WORD PAIRS GROUPS 
        [sameWordGroup, reverseWordGroup, probeWordGroup, targetWordGroup, diffWordGroup] = createAcrossWordGroups(firstwordpairs, secondwordpairs);
        
        % sessionblock directories for block(i) and block(i+1)
        sessionFirstBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        sessionSecondBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock+1});
        
        %%- BUILD FEATURE MAT events X features X time
        [samePairFeatureMat1, samePairFeatureMat2] = buildAcrossDiffPairFeatureMat(sameWordGroup, sessionFirstBlockDir, sessionSecondBlockDir);
        [reversePairFeatureMat1, reversePairFeatureMat2] = buildAcrossDiffPairFeatureMat(reverseWordGroup, sessionFirstBlockDir, sessionSecondBlockDir);
        [diffPairFeatureMat1, diffPairFeatureMat2] = buildAcrossDiffPairFeatureMat(diffWordGroup, sessionFirstBlockDir, sessionSecondBlockDir);
        [probePairFeatureMat1, probePairFeatureMat2] = buildAcrossDiffPairFeatureMat(probeWordGroup, sessionFirstBlockDir, sessionSecondBlockDir);
        [targetPairFeatureMat1, targetPairFeatureMat2] = buildAcrossDiffPairFeatureMat(targetWordGroup, sessionFirstBlockDir, sessionSecondBlockDir);
        
        % reshape matrix to events X time X features
        samePairFeatureMat1 = permute(samePairFeatureMat1, [1 3 2]);
        samePairFeatureMat2 = permute(samePairFeatureMat2, [1 3 2]);
        reversePairFeatureMat1 = permute(reversePairFeatureMat1, [1 3 2]);
        reversePairFeatureMat2 = permute(reversePairFeatureMat2, [1 3 2]);
        diffPairFeatureMat1 = permute(diffPairFeatureMat1, [1 3 2]);
        diffPairFeatureMat2 = permute(diffPairFeatureMat2, [1 3 2]);
        probePairFeatureMat1 = permute(probePairFeatureMat1, [1 3 2]);
        probePairFeatureMat2 = permute(probePairFeatureMat2, [1 3 2]);
        targetPairFeatureMat1 = permute(targetPairFeatureMat1, [1 3 2]);
        targetPairFeatureMat2 = permute(targetPairFeatureMat2, [1 3 2]);
        
        minSampleSize = min([size(reversePairFeatureMat1,1), size(samePairFeatureMat1,1), ...
                            size(targetPairFeatureMat1,1), size( probePairFeatureMat1, 1), ...
                            size(diffPairFeatureMat1,1)]);
                        
        randIndices = randsample(size(reversePairFeatureMat1, 1), minSampleSize);
        reversePairFeatureMat1 = reversePairFeatureMat1(randIndices, :, :);
        reversePairFeatureMat2 = reversePairFeatureMat2(randIndices, :, :);
        randIndices = randsample(size(diffPairFeatureMat1, 1), minSampleSize);
        diffPairFeatureMat1 = diffPairFeatureMat1(randIndices, :, :);
        diffPairFeatureMat2 = diffPairFeatureMat2(randIndices, :, :);
        randIndices = randsample(size(samePairFeatureMat1, 1), minSampleSize);
        samePairFeatureMat1 = samePairFeatureMat1(randIndices, :, :);
        samePairFeatureMat2 = samePairFeatureMat2(randIndices, :, :);
        
        randIndices = randsample(size(targetPairFeatureMat1, 1), minSampleSize);
        targetPairFeatureMat1 = targetPairFeatureMat1(randIndices, :, :);
        targetPairFeatureMat2 = targetPairFeatureMat2(randIndices, :, :);
        randIndices = randsample(size(probePairFeatureMat1, 1), minSampleSize);
        probePairFeatureMat1 = probePairFeatureMat1(randIndices, :, :);
        probePairFeatureMat2 = probePairFeatureMat2(randIndices, :, :);
        
        size(samePairFeatureMat1)
        size(reversePairFeatureMat1)
        size(diffPairFeatureMat1)
        
        %%- BUILD REINSTATEMENT MATRICES
        % same Pairs
        [eventSame, featureSame] = compute_reinstatement(samePairFeatureMat1, samePairFeatureMat2);
        [eventDiff, featureDiff] = compute_reinstatement(diffPairFeatureMat1, diffPairFeatureMat2);
        [eventReverse, featureReverse] = compute_reinstatement(reversePairFeatureMat1, reversePairFeatureMat2);
        [eventProbe, featureProbe] = compute_reinstatement(probePairFeatureMat1, probePairFeatureMat2);
        [eventTarget, featureTarget] = compute_reinstatement(targetPairFeatureMat1, targetPairFeatureMat2);
        
        size(squeeze(mean(eventDiff(:, :, :),1)))
        size(eventDiff)
        size(featureTarget)
        size(featureDiff)
        
        % rand sample down the different word pair feature mat -> match
        % size
%         minSampleSize = min([size(eventSame,1), size(eventDiff,1), ...
%                             size(eventReverse,1), size(eventProbe, 1), ...
%                             size(eventTarget,1)]);
%         
%         randIndices = randsample(size(eventSame,1), minSampleSize);
%         eventSame = eventSame(randIndices,:,:);
%         randIndices = randsample(size(eventDiff,1), minSampleSize);
%         eventDiff = eventDiff(randIndices,:,:);
%         randIndices = randsample(size(eventReverse,1), minSampleSize);
%         eventReverse = eventReverse(randIndices,:,:);
%         randIndices = randsample(size(eventProbe,1), minSampleSize);
%         eventProbe = eventProbe(randIndices,:,:);
%         randIndices = randsample(size(eventTarget,1), minSampleSize);
%         eventTarget = eventTarget(randIndices,:,:);
        
        %%- Save Mat files
        figureFile = strcat(figureDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs',num2str(blocks{iBlock+1}));
        matFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs',num2str(blocks{iBlock+1}));
        %%- save reinstatement matrices
        save(strcat(matFile, '.mat'), 'eventSame', 'featureSame', ...
                                        'eventReverse', 'featureReverse', ...
                                        'eventDiff', 'featureDiff', ...
                                        'eventTarget', 'featureTarget', ...
                                        'eventProbe', 'featureProbe');
        
        % write debugging output to .txt file
        sameWordGroup = [sameWordGroup{:}];
        reverseWordGroup = [reverseWordGroup{:}];
        diffWordGroup = [diffWordGroup{:}];
        probeWordGroup = [probeWordGroup{:}];
        targetWordGroup = [targetWordGroup{:}];
        logFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs',num2str(blocks{iBlock+1}), '.txt');
        fid = fopen(logFile, 'w');
        fprintf(fid, '%6s \n', 'Block(i) word pairs:');
        fprintf(fid, '%6s \n', firstwordpairs{:});
        fprintf(fid, '\n %s \n', 'Block(i+1) word pairs:');
        fprintf(fid, '%6s \n', secondwordpairs{:});
        fprintf(fid, '\n %s \n', 'Same Pair Group:');
        fprintf(fid, '%6s vs %6s \n', sameWordGroup{:});
        fprintf(fid, '\n %s \n', 'Reverse Pair Group:');
        fprintf(fid, '%6s vs %6s \n', reverseWordGroup{:});
        fprintf(fid, '\n %s \n', 'Different Pair Group:');
        fprintf(fid, '%6s vs %6s \n', diffWordGroup{:});
        fprintf(fid, '\n %s \n', 'Probe Overlap Pair Group:');
        fprintf(fid, '%6s vs %6s \n', probeWordGroup{:});
        fprintf(fid, '\n %s \n', 'Target Overlap Pair Group:');
        fprintf(fid, '%6s vs %6s \n', targetWordGroup{:});
        
        %%- Plotting
        fig = {};
        figure
        fig{end+1} = subplot(321)
        imagesc(squeeze(mean(eventSame(:, :, :),1)));
        title(['Same Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(eventDiff,1)), ' events'])
        colorbar();
        clim = get(gca, 'clim');
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        
        fig{end+1} = subplot(323);
        imagesc(squeeze(mean(eventDiff(:, :, :),1)));
        title(['Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1), ...
            ' vs ', num2str(iBlock)])
        colorbar();
        set(gca, 'clim', clim);
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        
        fig{end+1} = subplot(325);
        imagesc(squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(eventDiff(:, :, :),1)));
        title({'Within-Blocks', ['Same-Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1), ...
            ' vs ', num2str(iBlock)]})
        colorbar();
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        
        fig{end+1} = subplot(322)
        imagesc(squeeze(mean(eventReverse(:, :, :),1)));
        title(['Reverse Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(eventReverse,1)), ' events'])
        colorbar();
        set(gca, 'clim', clim);
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        
        fig{end+1} = subplot(324);
        imagesc(squeeze(mean(eventProbe(:, :, :),1)));
        title(['Probe Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(eventProbe,1)), ' events'])
        colorbar();
        set(gca, 'clim', clim);
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        
        fig{end+1} = subplot(326);
        imagesc(squeeze(mean(eventTarget(:, :, :),1)));
        title(['Target Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(eventTarget,1)), ' events'])
        colorbar();
        set(gca, 'clim', clim);
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        %%- SAVE THE FIGURE AFTER CHANGING IT
        fig = gcf;
        fig.PaperUnits = 'inches';
        pos = [0.35, 3.65, 12.55, 7.50];
        fig.PaperPosition = pos;
        
        %%- Save Image
        print(figureFile, '-dpng', '-r0')
        savefig(figureFile)
        
        pause(0.1);
    end
end
end
