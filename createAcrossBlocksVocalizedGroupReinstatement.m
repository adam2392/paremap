function createAcrossBlocksVocalizedGroupReinstatement(subj, typeTransform, timeLock, referenceType)
    close all;
    clear all;
    clc;
    
%% PARAMETERS FOR RUNNING PREPROCESS
% subj = 'NIH034';
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
%%------------------ STEP 2: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(dataDir, TYPE_TRANSFORM, 'vocalization_sessiontargetwords')
sessions = dir(dataDir);
sessions = {sessions(3:end).name};

if strcmp(subj, 'NIH039')
    sessions = sessions([1,2,4]);
elseif strcmp(subj, 'NIH034')
    sessions = sessions([3, 4]);
end
sessions

blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};
    
allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                    'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                    'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', 'GLASS_JUICE', ...
                    'GLASS_GLASS', 'JUICE_JUICE'};

% saving figures dir.
if strcmp(typeTransform, 'morlet')
    ticks = [6:10:46];
    labels = [-2:1:3];
    timeZero = 26;
elseif strcmp(typeTransform, 'multitaper')
    ticks = [1:2:11];
    labels = [-2:1:3];
    timeZero = 5;
end
figureDir = strcat('./Figures/', subj, '/reinstatement/', TYPE_TRANSFORM,'/across_blocks_vocalizationWord/');
matDir = strcat('./Figures/', subj, '/reinstatement_mat/', TYPE_TRANSFORM,'/across_blocks_vocalizationWord/');
LT = 1.5 %line thickness

if ~exist(figureDir, 'dir')
    mkdir(figureDir)
end
if ~exist(matDir, 'dir')
    mkdir(matDir)
end

%% RUN ANALYSIS
%%- LOOP THROUGH SESSIONS AND BLOCKS
for iSesh=1:length(sessions),
    for iBlock=1:length(blocks)-1, % loop through first 5 blocks and do across blocks analysis
        % initialize feature matrix cell
        allVocalizedIndices = zeros(length(allVocalizedPairs), 1);

        %%- 01: BUILD WORD PAIRS
        % get word pairs in this session-block
        targetWordsOne = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        targetWordsOne = {targetWordsOne(3:end).name};
        targetWordsTwo = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock+1}));
        targetWordsTwo = {targetWordsTwo(3:end).name};
        wordPairs = createAcrossVocalizedWordGroups(targetWordsOne, targetWordsTwo);

        %%- BUILD FEATURE MATRIX
        %%- set plotting directories and meta data output
        figureFile = strcat(figureDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs', num2str(blocks{iBlock+1}));
        matFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs', num2str(blocks{iBlock+1}));  

        % meta data logging
        logFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs', num2str(blocks{iBlock+1}), '.txt');
        fid = fopen(logFile, 'w');
        fprintf(fid, '%6s \n', 'Block(i) Target Words (vocalized):');
        fprintf(fid, '%6s \n', targetWordsOne{:}); 
        fprintf(fid, '%6s \n', 'Block(i+1) Target Words (vocalized):');
        fprintf(fid, '%6s \n', targetWordsTwo{:}); 
        fprintf(fid, '%6s \n', 'Word Pairs:');
        fprintf(fid, '%6s \n', wordPairs{:}); 

        sessionFirstBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        sessionSecondBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock+1});

        eventReinMat = {};
        featureReinMat = {};
        %%- for a certain session, block's vocalized words (targetWords)
        % LOOP through all wordPairs and build feature mat
        for iPair=1:length(wordPairs)
            wordSplit = strsplit(wordPairs{iPair}, '_');
            firstWord = wordSplit{1};
            secondWord = wordSplit{2};
            
            %%- 02: BUILD FEATURE MATRIX 
            [pairFeatureMat1, pairFeatureMat2] = buildAcrossPairFeatureMat(firstWord, secondWord, sessionFirstBlockDir, sessionSecondBlockDir);                
            % change time/features dimensions
            pairFeatureMat1 = permute(pairFeatureMat1, [1 3 2]);
            pairFeatureMat2 = permute(pairFeatureMat2, [1 3 2]);

            %%- 03: BUILD REINSTATEMENT MATRICES
            [eventRein, featureRein] = compute_reinstatement(pairFeatureMat1, pairFeatureMat2);
            size(eventRein)
            size(featureRein)

            %%- 02: DETERMINE INDEX IN OUR LIST OF WORD PAIRS
            allVocalizedPairs;
            checkOne = strjoin({firstWord, secondWord}, '_');
            checkTwo = strjoin({secondWord, firstWord}, '_');

            %%- FIND INDEX IN MASTER LIST OF WORDPAIRS
            % Check if this pair is in our list of 15 vocalization pairs
            if (ismember(checkOne, allVocalizedPairs) ||...
                ismember(checkTwo, allVocalizedPairs))

                %- find first combintation
                index = cellfun(@(x) strcmp(checkOne, x), allVocalizedPairs, 'UniformOutput', 0);
                if isempty(find([index{:}] == 1)) %- find second combination
                    index = cellfun(@(x) strcmp(checkTwo, x), allVocalizedPairs, 'UniformOutput', 0);
                end
                index = find([index{:}] == 1);

                %- neither combination, then it is an incorrec vocalized
                %word
                if ~strcmp(allVocalizedPairs{index}, checkOne) &&...
                        ~strcmp(allVocalizedPairs{index}, checkTwo)
                    disp('error?');
                end
            end

            %%- 03: BUILD FEATURE MATRIX UP USING CELL ARRAY
            % build onto feature matrices in cell mat
            if allVocalizedIndices(index) == 0
                eventReinMat{index} = eventRein;
                featureReinMat{index} = featureRein;
            else
                eventReinMat{index} = cat(1, eventReinMat{index}, eventRein);
                featureReinMat{index} = cat(1, featureReinMat{index}, featureRein);
            end

            allVocalizedIndices(index) = 1;
        end %loop through word pairs -> built feature matrix

        %%- SAVE MAT FILE PER BLOCK
        save(strcat(matFile, '.mat'), 'eventReinMat', 'featureReinMat', 'allVocalizedIndices');
        
        %%- 04: PLOTTING
        fig = figure;
        clim = [0 0]; %initialize colorbar
        fa = {};
        
        for iPlot=1:length(eventReinMat)
            iPlot
            eventRein = eventReinMat{iPlot};
            wordSplit = strsplit(allVocalizedPairs{iPlot}, '_');
            wordone = wordSplit{1};
            wordtwo = wordSplit{2};

            if allVocalizedIndices(iPlot) == 1
                fa{iPlot} = subplot(4, 4, iPlot);
                imagesc(squeeze(mean(eventRein(:,:,:),1)));
                title({['Cosine Similarity for Block ', num2str(iBlock-1), ...
                    'vs', num2str(iBlock), ' for '], [wordone, ' vs ', wordtwo, ...
                    ' (', num2str(size(eventRein,1)), ' events)']});
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
                colorbar();
                tempclim = get(gca, 'clim');
                clim(1) = min(tempclim(1), clim(1));
                clim(2) = max(tempclim(2), clim(2));
                plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
                plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
            else
                disp('did nothing')
            end
        end
        % change the color limit to the max in the group for comparison
        for i=1:length(allVocalizedPairs)
            if allVocalizedIndices(i) == 1
                fa{i}.CLim = clim;
            end
        end

        % change figure dimensions before saving
        fig = gcf;
        fig.PaperUnits = 'inches';
        pos = [5.1667 0.6806 9.9722 10.3889];
        fig.PaperPosition = pos;

        %%- Save the image
        print(figureFile, '-dpng', '-r0')
        savefig(figureFile)

        close all
    end %loop thru blocks
end %loop thru sessions
% end