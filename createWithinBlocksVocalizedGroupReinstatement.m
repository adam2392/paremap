function createWithinBlocksVocalizedGroupReinstatement(subj, typeTransform, timeLock, referenceType)    
clc;
close all;
    
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
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

% all target word comparisons that exist 
allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                    'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                    'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', 'GLASS_JUICE', ...
                    'GLASS_GLASS', 'JUICE_JUICE'};
    
%%- SAVING FIGURES OPTIONS
figureDir = strcat('./Figures/', subj, '/reinstatement/', TYPE_TRANSFORM, '/within_blocks_vocalizationWord/');
matDir = strcat('./Figures/', subj, '/reinstatement_mat/', TYPE_TRANSFORM, '/within_blocks_vocalizationWord/');         
if ~exist(figureDir, 'dir') mkdir(figureDir); end
if ~exist(matDir, 'dir')    mkdir(matDir);    end

LT = 1.5; % line thickness

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
            
%% CREATE VOCALIZED WORD GROUPS
%%- LOOP THROUGH SESSIONS AND BLOCKS
%%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:length(blocks),
        fprintf('%6s \n', strcat('On session ', num2str(iSesh), ' and block ', num2str(iBlock)));
        
        % get word pairs in this session-block
        targetWords = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        targetWords = {targetWords(3:end).name};
        % get all word pairs consisting of the target words
        wordPairs = createWithinVocalizedWordGroups(targetWords);
        
        %%- set plotting directories and meta data output
        figureFile = strcat(figureDir, sessions{iSesh}, '_', blocks{iBlock});
        matFile = strcat(matDir, sessions{iSesh}, '_', blocks{iBlock});
        %%- do some error logging in a .txt file
        wordpairs = [wordPairs{:}];
        logFile = strcat(matDir, sessions{iSesh}, '_', blocks{iBlock}, '.txt');
        fid = fopen(logFile, 'w');
        fprintf(fid, '%6s \n', 'Block(i) Target Words (vocalized):');
        fprintf(fid, '%6s \n', fullfile(dataDir, sessions{iSesh}, blocks{iBlock})); 
        fprintf(fid, '%6s \n', 'Word Pairs:');
        fprintf(fid, '%6s \n', wordpairs{:}); 

        % initialize cell arrays to store data
        eventReinMat = {};
        featureReinMat = {};
        sessionBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        
        % a vector of flags to determine whether or not a comparison was
        % already made between two vocalized words
        allVocalizedIndices = zeros(length(allVocalizedPairs), 1);

        %%- NOW LOOP THROUGH ALL PAIRS OF TARGET WORDS TO COMPARE REINSTATEMENT MAPS
        for iWord=1:length(wordPairs),            
            wordone = wordPairs{iWord}{1}; % first vocalized word
            wordtwo = wordPairs{iWord}{2}; % second vocalized word

            %%- 02: BUILD FEATURE MATRICES
            %%- get pair feature matrices for every vocalized word comparison
            [pairFeatureMat1, pairFeatureMat2] = buildWithinPairFeatureMat(wordone, wordtwo, sessionBlockDir);
            % change time/features dimensions -> events X time X features
            pairFeatureMat1 = permute(pairFeatureMat1, [1 3 2]);
            pairFeatureMat2 = permute(pairFeatureMat2, [1 3 2]);

            %%- 03: BUILD REINSTATEMENT MATRICES
            [eventRein, featureRein] = compute_reinstatement(pairFeatureMat1, pairFeatureMat2);
%             size(pairFeatureMat1)
%             size(pairFeatureMat2)
%             size(eventRein)
%             size(featureRein)
            checkOne = strjoin({wordone, wordtwo}, '_');
            checkTwo = strjoin({wordtwo, wordone}, '_');
            disp('Finished building reinstatement matrices for: ');
            checkOne
            checkTwo
            %%- 04: FIND INDEX IN MASTER LIST OF WORDPAIRS
            % Check if this pair is in our list of 15 vocalization pairs
            if (ismember(checkOne, allVocalizedPairs) ||...
                ismember(checkTwo, allVocalizedPairs))

                %- find first combintation
                ind = cellfun(@(x) strcmp(checkOne, x), allVocalizedPairs, 'UniformOutput', 0);
                if isempty(find([ind{:}] == 1)) %- find second combination
                    ind = cellfun(@(x) strcmp(checkTwo, x), allVocalizedPairs, 'UniformOutput', 0);
                end
                ind = find([ind{:}] == 1);

                %- neither combination, then it is an incorrec vocalized
                %word
                if ~strcmp(allVocalizedPairs{ind}, checkOne) &&...
                        ~strcmp(allVocalizedPairs{ind}, checkTwo)
                    disp('error?');
                end
            end

            %%- 05: BUILD FEATURE MATRIX PER VOCALIZED WORD PAIR
            if allVocalizedIndices(ind) == 0
                eventReinMat{ind} = eventRein;
                featureReinMat{ind} = featureRein;
            else
                eventReinMat{ind} = cat(1, eventReinMat{ind}, eventRein);
                featureReinMat{ind} = cat(1, featureReinMat{ind}, featureRein);
            end

            allVocalizedIndices(ind) = 1;
        end
        size(eventReinMat)

        %%- finished building reinstatement map of within-block targetWord
        %%Pairs, and now save per session/block
        %%- SAVE MAT FILE PER BLOCK
        save(strcat(matFile, '.mat'), 'eventReinMat', 'featureReinMat', 'allVocalizedIndices');

        %%- 04: PLOTTING
        fig = figure;
        clim = [0 0]; %initialize colorbar
        fa = {};
        for iPlot=1:length(eventReinMat)
            eventRein = eventReinMat{iPlot};
            wordSplit = strsplit(allVocalizedPairs{iPlot}, '_');
            wordone = wordSplit{1};
            wordtwo = wordSplit{2};

            if allVocalizedIndices(iPlot) == 1
                fa{iPlot} = subplot(4, 4, iPlot);
                imagesc(squeeze(mean(eventRein(:,:,:),1)));
                title({['Cosine Similarity for Block ', num2str(iBlock-1), ...
                    ' for '], [wordone, ' vs ', wordtwo, ' (', num2str(size(eventRein,1)), ' events)']});
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
        end % end of loop through words            
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
    end % loop through blocks
end % loop through sessions

