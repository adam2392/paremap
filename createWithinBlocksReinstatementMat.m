%         -- Creating matrices that are eventsXtimeXfeatures, which can be
%         fed into compute_reinstatement.m
%         -- This is done for within blocks analysis of the paremap task
%        

function createWithinBlocksReinstatementMat(subj, typeTransform, timeLock, referenceType)
close all;
clc;

%% PARAMETERS FOR RUNNING PREPROCESS
% subj = 'NIH039';
% timeLock = 'vocalization';
% referenceType = 'bipolar';
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
%%------------------ STEP 1: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(dataDir, TYPE_TRANSFORM, CUE_LOCK)
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
sessions

blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

%%- SAVING FIGURES OPTIONS
if strcmp(CUE_LOCK, 'vocalization')
    figureDir = strcat('./Figures/', subj, '/reinstatement/', TYPE_TRANSFORM,'/within_blocks_vocalization/');
    matDir = strcat('./Figures/', subj, '/reinstatement_mat/', TYPE_TRANSFORM,'/within_blocks_vocalization/');
elseif strcmp(CUE_LOCK, 'matchword')
    figureDir = strcat('./Figures/', subj, '/reinstatement/', TYPE_TRANSFORM,'/within_blocks_matchword/');
    matDir = strcat('./Figures/', subj, '/reinstatement_mat/', TYPE_TRANSFORM,'/within_blocks_matchword/');
elseif strcmp(CUE_LOCK, 'probeword')
    figureDir = strcat('./Figures/', subj, '/reinstatement/', TYPE_TRANSFORM,'/within_blocks_probeon/');
    matDir = strcat('./Figures/', subj, '/reinstatement_mat/', TYPE_TRANSFORM,'/within_blocks_probeon/');
end

% make directory, if doens't already exist
if ~exist(figureDir) mkdir(figureDir); end
if ~exist(matDir)    mkdir(matDir);    end
% set linethickness
LT = 1.5;

% load in an example file to get the -> labels, ticks and timeZero
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
    for iBlock=1:length(blocks),
        fprintf('%6s \n', strcat('On session ', num2str(iSesh), ' and block ', num2str(iBlock)));
        
        % get word pairs in this session-block
        wordpairs = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        wordpairs = {wordpairs(3:end).name};
        
        %%- CREATE WITHIN WORD PAIRS GROUPS 
        [sameWordGroup, reverseWordGroup, diffWordGroup] = createWithinWordGroups(wordpairs);
        
        %%- BUILD FEATURE MATRIX 
        fprintf('Building pair feature matrices for same,reverse,diff wordGroup \n');
        sessionBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        [samePairFeatureMat1, samePairFeatureMat2] = buildSamePairFeatureMat(sameWordGroup, sessionBlockDir);
        [reversePairFeatureMat1, reversePairFeatureMat2] = buildDiffPairFeatureMat(reverseWordGroup, sessionBlockDir);
        [diffPairFeatureMat1, diffPairFeatureMat2] = buildDiffPairFeatureMat(diffWordGroup, sessionBlockDir);
        
        % reshape matrix to events X time X features
        samePairFeatureMat1 = permute(samePairFeatureMat1, [1 3 2]);
        samePairFeatureMat2 = permute(samePairFeatureMat2, [1 3 2]);
        reversePairFeatureMat1 = permute(reversePairFeatureMat1, [1 3 2]);
        reversePairFeatureMat2 = permute(reversePairFeatureMat2, [1 3 2]);
        diffPairFeatureMat1 = permute(diffPairFeatureMat1, [1 3 2]);
        diffPairFeatureMat2 = permute(diffPairFeatureMat2, [1 3 2]);
        
        randIndices = randsample(size(reversePairFeatureMat1, 1), size(samePairFeatureMat1, 1));
        reversePairFeatureMat1 = reversePairFeatureMat1(randIndices, :, :);
        reversePairFeatureMat2 = reversePairFeatureMat2(randIndices, :, :);
        randIndices = randsample(size(diffPairFeatureMat1, 1), size(samePairFeatureMat1, 1));
        diffPairFeatureMat1 = diffPairFeatureMat1(randIndices, :, :);
        diffPairFeatureMat2 = diffPairFeatureMat2(randIndices, :, :);
        
        size(samePairFeatureMat1)
        size(reversePairFeatureMat1)
        size(diffPairFeatureMat1)
        
        %%- BUILD REINSTATEMENT MATRICES
        % same Pairs
        fprintf('Computing reinstatement matrices \n');
        [eventSame, featureSame] = compute_reinstatement(samePairFeatureMat1, samePairFeatureMat2);
        [eventReverse, featureReverse] = compute_reinstatement(reversePairFeatureMat1, reversePairFeatureMat2);
        [eventDiff, featureDiff] = compute_reinstatement(diffPairFeatureMat1, diffPairFeatureMat2);
        
        size(squeeze(mean(eventDiff(:, :, :),1)))
        size(eventSame)
        size(featureDiff)
        size(featureSame)
        
%         size(samePairFeatureMat1)
%         size(samePairFeatureMat2)
%         size(diffPairFeatureMat1)
%         size(diffPairFeatureMat2)

        % rand sample down the different word pair feature mat -> match
        % size
%         randIndices = randsample(size(eventDiff,1), size(eventSame,1));
%         eventDiff = eventDiff(randIndices,:,:);
%         randIndices = randsample(size(eventReverse,1), size(eventSame,1));
%         eventReverse = eventReverse(randIndices,:,:);
        
        %%- Save Mat Files
        figureFile = strcat(figureDir, sessions{iSesh}, '-', num2str(blocks{iBlock}));
        matFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}));
        %%- save reinstatement matrices
        save(strcat(matFile, '.mat'), 'eventSame', 'featureSame', ...
                                        'eventReverse', 'featureReverse',...
                                        'eventDiff', 'featureDiff');
        
        % write debugging output to .txt file
        sameWordGroup = [sameWordGroup{:}];
        reverseWordGroup = [reverseWordGroup{:}];
        diffWordGroup = [diffWordGroup{:}];
        logFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), '.txt');
        fid = fopen(logFile, 'w');
        fprintf(fid, '%6s \n', 'Block(i) word pairs:');
        fprintf(fid, '%6s \n', wordpairs{:});
        fprintf(fid, '\n %s \n', 'Same Pair Group:');
        fprintf(fid, '%6s vs %6s \n', sameWordGroup{:});
        fprintf(fid, '\n %s \n', 'Reverse Pair Group:');
        fprintf(fid, '%6s vs %6s \n', reverseWordGroup{:});
        fprintf(fid, '\n %s \n', 'Different Pair Group:');
        fprintf(fid, '%6s vs %6s \n', diffWordGroup{:});                          
                                                     
        %%- Plotting
        fig = {};
        figure;
        fig{end+1} = subplot(221);
        imagesc(squeeze(mean(eventSame(:, :, :),1)));
        title(['Same Pairs Cosine Similarity for Block ', num2str(iBlock-1) ...
            ' with ', num2str(size(eventDiff,1)), ' events'])
        colorbar();
        clim = get(gca, 'clim');
        axis square
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)

        fig{end+1} = subplot(222);
        imagesc(squeeze(mean(eventDiff(:, :, :),1)));
        title(['Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
        colorbar();
        set(gca, 'clim', clim);
        axis square
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)

        fig{end+1} = subplot(223);
        imagesc(squeeze(mean(eventReverse(:, :, :),1)));
        title(['Reverse Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
        colorbar();
        set(gca, 'clim', clim);
        axis square
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)

        fig{end+1} = subplot(224);
        imagesc(squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(eventDiff(:, :, :),1)));
        title({'Within-Blocks', ['Same-Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)]})
        colorbar();
        axis square
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)

        %%- Save Image
        print(figureFile, '-dpng', '-r0')
        savefig(figureFile)
        
        pause(0.1); 
    end % loop through blocks
end % loop through sessions
end
