% function createAcrossBlocksVocalizedGroupReinstatement(subj)
    close all;
    clear all;
    clc;
    
    subj = 'NIH039';
    sessNum = [0, 1, 2];
    addpath('./m_reinstatement/');
    
    %% INITIALIZATION OF SESSION AND BLOCKS TO LOOK AT
    TYPE_TRANSF = 'morlet_spec_vocalization';
    disp('WITHIN BLOCKS');
    disp(TYPE_TRANSF);
    
    dataDir = strcat('./condensed_data_', subj);
    dataDir = fullfile(dataDir, TYPE_TRANSF);
    sessions = dir(dataDir);
    sessions = {sessions(3:end).name};
    % sessions = sessions(3:end);

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

    %% RUN ANALYSIS
    %%- LOOP THROUGH SESSIONS AND BLOCKS
    for iSesh=1:length(sessions),
        for iBlock=1:length(blocks)-1, % loop through first 5 blocks and do across blocks analysis
            % initialize feature matrix cell
%             eventReinMat = cell(length(allVocalizedPairs), 1);
%             featureReinMat = cell(length(allVocalizedPairs), 1);
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
            figureDir = strcat('./Figures/', subj, '/reinstatement/across_blocks_vocalizationWord/');
            matDir = strcat('./Figures/', subj, '/reinstatement_mat/across_blocks_vocalizationWord/');
            figureFile = strcat(figureDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs', num2str(blocks{iBlock+1}));
            matFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs', num2str(blocks{iBlock+1}));  
            if ~exist(figureDir)
                mkdir(figureDir)
            end
            if ~exist(matDir)
                mkdir(matDir)
            end
            
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
            
            % LOOP through all wordPairs and build feature mat
            for iPair=1:length(wordPairs)
                wordSplit = strsplit(wordPairs{iPair}, '_');
                firstWord = wordSplit{1};
                secondWord = wordSplit{2};
                
                wordSplit
                %%- BUILD FEATURE MATRIX 
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

                    index = cellfun(@(x) strcmp(checkOne, x), allVocalizedPairs, 'UniformOutput', 0);
                    if isempty(find([index{:}] == 1))
                        index = cellfun(@(x) strcmp(checkTwo, x), allVocalizedPairs, 'UniformOutput', 0);
                    end
                    index = find([index{:}] == 1);

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
                
                % save the relevant mat files
                if ~exist(strcat(matFile, '.mat'))
                    fieldname = strcat(firstWord, '_', secondWord);
                    eventReinData.(fieldname) = eventRein;
                    featureReinData.(fieldname) = featureRein;
                else
                    load(strcat(matFile, '.mat'));
                    fieldname = strcat(firstWord, '_', secondWord);
                    eventReinData.(fieldname) = eventRein;
                    featureReinData.(fieldname) = featureRein;
                end
                save(strcat(matFile, '.mat'), 'eventReinData', 'featureReinData');
                
                allVocalizedIndices(index) = 1;
            end %loop through word pairs -> built feature matrix
    
            %%- 04: PLOTTING
            fig = figure;
            clim = [0 0]; %initialize colorbar
            fa = {};
            LT = 1.5 %line thickness
            
            timeZero = 16; % for 0 seconds from -2->4
            ticks = [6:10:56]; % for 6 seconds of data
            labels = [-1:1:4];
            
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