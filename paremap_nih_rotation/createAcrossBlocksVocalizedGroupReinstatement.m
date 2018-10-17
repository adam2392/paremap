function createAcrossBlocksVocalizedGroupReinstatement(subj, typeTransform, referenceType)
close all;
%     clear all;
clc;

% subj = 'NIH034';
% timeLock = 'vocalization';
% referenceType = 'bipolar';
% typeTransform = 'morlet';


addpath('./reinstatement_vocalizedWords/');
%% PARAMETERS FOR RUNNING PREPROCESS
expected_timeLocks = {'vocalization', 'matchword', 'probeword'};
expected_transforms = {'morlet', 'multitaper'};
REF_TYPES = {'noreref', 'bipolar', 'global'};
if ~ismember(referenceType, REF_TYPES)
    disp('reference types are noreref, bipolar, or global');
end
if ~ismember(typeTransform, expected_transforms)
    disp('transform types are morlet, multitaper');
end
THIS_REF_TYPE = referenceType; 
TYPE_TRANSFORM = strcat(typeTransform, '_', referenceType);

addpath('./m_reinstatement/');
    
%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjDataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(subjDataDir, strcat(typeTransform, '_', referenceType, '_', 'vocalization_sessiontargetwords'));
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};
sessions

% all target word comparisons that exist. A_B, vocalizations of A is compared with vocalizations of B
allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                    'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                    'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', ...
                    'GLASS_JUICE', 'GLASS_GLASS', ...
                    'JUICE_JUICE'};

% saving figures dir.
figureDir = strcat('./Figures/', subj, '/reinstatement/',typeTransform, '_', referenceType, '_', 'across_blocks_vocalizationWord/');
matDir = strcat('./Figures/', subj, '/reinstatement_mat/', typeTransform, '_', referenceType, '_','across_blocks_vocalizationWord/');
if ~exist(figureDir, 'dir'), mkdir(figureDir); end
if ~exist(matDir, 'dir'), mkdir(matDir); end

LT = 1.5 %line thickness
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

%% RUN ANALYSIS
%%- LOOP THROUGH SESSIONS AND BLOCKS
for iSesh=1:length(sessions),
    for iBlock=1:length(blocks)-1, % loop through first 5 blocks and do across blocks analysis
        fprintf('%6s \n', strcat('On session ', num2str(iSesh), ' and block ', num2str(iBlock)));
        
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

        % session-block directories for each targetWord
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

            %%- 02: DETERMINE INDEX IN OUR LIST OF WORD PAIRS
            checkOne = strjoin({firstWord, secondWord}, '_');
            checkTwo = strjoin({secondWord, firstWord}, '_');
%             checkOne
%             checkTwo

            %%- FIND INDEX IN MASTER LIST OF WORDPAIRS
            % Check if this pair is in our list of 15 vocalization pairs
            if (ismember(checkOne, allVocalizedPairs) ||...
                ismember(checkTwo, allVocalizedPairs))

                %- find first combintation
                ind = cellfun(@(x) strcmp(checkOne, x), allVocalizedPairs, 'UniformOutput', 0);
                if isempty(find([ind{:}] == 1)) %- find second combination
                    ind = cellfun(@(x) strcmp(checkTwo, x), allVocalizedPairs, 'UniformOutput', 0);
                end
                ind = find([ind{:}] == 1);

                %- neither combination, then it is an incorrect vocalized
                %word
                if ~strcmp(allVocalizedPairs{ind}, checkOne) &&...
                        ~strcmp(allVocalizedPairs{ind}, checkTwo)
                    disp('error?');
                end
            end
%             ind
            %%- 03: BUILD FEATURE MATRIX UP USING CELL ARRAY
            % build onto feature matrices in cell mat
            if allVocalizedIndices(ind) == 0
                eventReinMat{ind} = eventRein;
                featureReinMat{ind} = featureRein;
            else
                eventReinMat{ind} = cat(1, eventReinMat{ind}, eventRein);
                featureReinMat{ind} = cat(1, featureReinMat{ind}, featureRein);
            end

            allVocalizedIndices(ind) = 1;
        end %loop through word pairs -> built feature matrix

        %%- SAVE MAT FILE PER BLOCK
        save(strcat(matFile, '.mat'), 'eventReinMat', 'featureReinMat', 'allVocalizedIndices');
        
        size(eventReinMat)
        
        %%- 04: PLOTTING
        fig = figure;
        clim = [0 0]; %initialize colorbar
        fa = {};
        
        %%- Loop and plot all targetWord comparisons
        for iPlot=1:length(eventReinMat)
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
        pos = [0    0.6667   17.5972   10.4028];
        fig.PaperPosition = pos;

        %%- Save the image
        print(figureFile, '-dpng', '-r0')
        savefig(figureFile)

        close all
    end %loop thru blocks
end %loop thru sessions
% end