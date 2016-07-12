function createWithinBlocksVocalizedGroupReinstatement(subj, typeTransform, timeLock, referenceType)    

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

%% LOAD EVENTS STRUCT AND SET DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirJhu = '/home/adamli/paremap';     % work
eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap'; 
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home

% Determine which directory we're working with automatically
if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

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
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                    'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                    'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', 'GLASS_JUICE', ...
                    'GLASS_GLASS', 'JUICE_JUICE'};
    
%%- SAVING FIGURES OPTIONS
if strcmp(typeTransform, 'morlet')
    ticks = [6:10:46];
    labels = [-2:1:3];
    timeZero = 26;
else
    timeZero = 5;
    ticks = [1:2:11];
    labels = [-2:1:3];
end
figureDir = strcat('./Figures/', subj, '/reinstatement/', TYPE_TRANSFORM, '/within_blocks_vocalizationWord/');
matDir = strcat('./Figures/', subj, '/reinstatement_mat/', TYPE_TRANSFORM, '/within_blocks_vocalizationWord/');         
LT = 1.5; % line thickness
if ~exist(figureDir, 'dir')
    mkdir(figureDir)
end
if ~exist(matDir, 'dir')
    mkdir(matDir)
end
                
%% CREATE VOCALIZED WORD GROUPS
%%- LOOP THROUGH SESSIONS AND BLOCKS
%%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:length(blocks),
        % get word pairs in this session-block
        targetWords = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        targetWords = {targetWords(3:end).name};
        % get all pairs of target words
        wordPairs = createWithinVocalizedWordGroups(targetWords);
        wordPairs{:}
        
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
        
        %%- NOW LOOP THROUGH ALL TARGETWORDS SAVED FOR SESSION/BLOCK
        allVocalizedIndices = zeros(length(allVocalizedPairs), 1);

        %%- NOW LOOP THROUGH ALL PAIRS OF TARGET WORDS TO COMPARE
        %%REINSTATEMENT MAPS
        for iWord=1:length(wordPairs),            
            wordone = wordPairs{iWord}{1}; % first vocalized word
            wordtwo = wordPairs{iWord}{2}; % second vocalized word

            %%- 02: BUILD FEATURE MATRICES
            %%- get pair feature matrices for every vocalized word comparison
            [pairFeatureMat1, pairFeatureMat2] = buildWithinPairFeatureMat(wordone, wordtwo, sessionBlockDir);
            % change time/features dimensions
            pairFeatureMat1 = permute(pairFeatureMat1, [1 3 2]);
            pairFeatureMat2 = permute(pairFeatureMat2, [1 3 2]);

            size(pairFeatureMat1)
            size(pairFeatureMat2)

            %%- 03: BUILD REINSTATEMENT MATRICES
            [eventRein, featureRein] = compute_reinstatement(pairFeatureMat1, pairFeatureMat2);
            size(eventRein)
            size(featureRein)


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
            iPlot
            eventRein = eventReinMat{iPlot};
            wordSplit = strsplit(allVocalizedPairs{iPlot}, '_');
            wordone = wordSplit{1};
            wordtwo = wordSplit{2};

            if allVocalizedIndices(iPlot) == 1
                fa{iWord} = subplot(5, 2, iWord);
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

