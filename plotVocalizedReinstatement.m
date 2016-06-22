close all
clc
clear all

%%- Section for Pulling reinstatement matrices already produced and just
%%make different plots
ANALYSIS_TYPE = {'within_blocks', 'across_blocks'};
EVENT_SYNC = {'probeon', 'vocalizationWord'};

subj = 'NIH039';
ANALYSIS = ANALYSIS_TYPE{2};
SYNC = EVENT_SYNC{2};

disp(ANALYSIS)
disp(SYNC)

% file dir for all the saved mat files
fileDir = strcat('./Figures/', subj, '/reinstatement_mat/', ANALYSIS, '_', SYNC, '/')


% load in an example data directory to get session names and block number
dataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(dataDir, 'morlet_spec');
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
if strcmp(subj, 'NIH039')
    sessions = sessions([1,2,4]);
elseif strcmp(subj, 'NIH034')
    sessions = sessions([3, 4]);
end
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

% set which blocks to analyze
if strcmp(ANALYSIS, 'across_blocks')
    lenBlocks = length(blocks)-1;
else
    lenBlocks = length(blocks);
end

sessions

%%- PLOTTING OPTIONS
% set linethickness
LT = 1.5;
if strcmp(SYNC, 'vocalization')
    ticks = [6:10:56];
    labels = [-3:1:2];
    timeZero = 36;
elseif strcmp(SYNC, 'probeon')
    ticks = [6:10:56];
    labels = [0:1:5];
    timeZero = 6;
end

allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                        'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                        'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', 'GLASS_JUICE', ...
                        'GLASS_GLASS', 'JUICE_JUICE'};

%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    % initialize average feature matrix cell
    eventReinMat = cell(length(allVocalizedPairs), 1);
    featureReinMat = cell(length(allVocalizedPairs), 1);
    allVocalizedIndices = zeros(length(allVocalizedPairs), 1);
    % Looking at the entire reinstatement
    figureFile = strcat('./Figures/', subj, '/reinstatement/', ANALYSIS, '_', SYNC, '/', ...
        'summary_', sessions{iSesh});
    
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:lenBlocks,
        if strcmp(ANALYSIS, 'across_blocks')
            sessionBlockName = strcat(sessions{iSesh}, '-', blocks{iBlock}, 'vs', blocks{iBlock+1});
        else
            sessionBlockName = strcat(sessions{iSesh}, '-', blocks{iBlock});
        end
        fileToLoad = fullfile(fileDir, sessionBlockName);
        data = load(fileToLoad);
        eventReinData = data.eventReinData;
        featureReinData = data.featureReinData;
        
        wordPairs = fieldnames(eventReinData);
        
        %%- Now loop through all these words and append to an averaged
        %%dictionary within cell array
        for iPair=1:length(wordPairs)
            %%- 01: Get the words in this session-block
            wordSplit = strsplit(wordPairs{iPair}, '_');
            firstWord = wordSplit{1};
            secondWord = wordSplit{2};
            checkOne = strjoin({firstWord, secondWord}, '_');
            checkTwo = strjoin({secondWord, firstWord}, '_');
            
            % extract the corresponding matrix from struct
            eventRein = eventReinData.(wordPairs{iPair});
            featureRein = featureReinData.(wordPairs{iPair});
            
            %%- 02: Find index in Master List
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
            
            %%- 03: BUILD AVERAGED FEATURE MATRIX using Cell Array
            if allVocalizedIndices(index) == 0
                eventReinMat{index} = eventRein;
                featureReinMat{index} = featureRein;
            else
                eventReinMat{index} = cat(1, eventReinMat{index}, eventRein);
                featureReinMat{index} = cat(1, featureReinMat{index}, featureRein);
            end
            
            allVocalizedIndices(index) = 1;
        end
    end % loop through blocks
    
    %%- 04: PLOTTING
    fig = figure;
    clim = [0 0]; %initialize colorbar
    fa = {};
    LT = 1.5 %line thickness

    timeZero = 16; % for 0 seconds from -2->4
    ticks = [6:10:56]; % for 6 seconds of data
    labels = [-1:1:4];
    
    %%- Plot each word pairing separately
    for iPlot=1:length(allVocalizedPairs)
        eventRein = eventReinMat{iPlot};
        
        % split words
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
    pos = [0    0.6806   20.0000   10.2917];
    fig.PaperPosition = pos;

    %%- Save the image
    print(figureFile, '-dpng', '-r0')
    savefig(figureFile)

    close all
end % loop through sessions