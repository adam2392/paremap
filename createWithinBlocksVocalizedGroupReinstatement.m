% function createWithinBlocksVocalizedGroupReinstatement(subj)
    close all;
    
    subj = 'NIH034';
    sessNum = [0, 1, 2];
    addpath('./m_reinstatement/');
    
    %% LOAD EVENTS STRUCT AND SET DIRECTORIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eegRootDirJhu = '/home/adamli/paremap';     % work
    eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap'; 
    eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home


    % Determine which directory we're working with automatically
    if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
    elseif length(dir(eegRootDirJhu))>0, eegRootDir = eegRootDirJhu;
    elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
    else   error('Neither Work nor Home EEG directories exist! Exiting'); end

    % Either go through all the sessions, or a specific session
    if sessNum == -1 | length(sessNum)>1, % all sessions
        disp('STEP 1: Going through all sessions')
        session = 'Meta Session [all]';
        behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap');
        sessStr = '[all]';
    else                                  % one session
        disp('STEP 1: Going through one session')
        session = sprintf('session_%d',sessNum);
        behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paRemap/', session);
        sessStr = sprintf('[%d]',sessNum);
    end

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
    
    length(allVocalizedPairs)
    
    %% CREATE VOCALIZED WORD GROUPS
    %%- LOOP THROUGH SESSIONS AND BLOCKS
    for iSesh=1:length(sessions),
        for iBlock=1:length(blocks),
            allVocalizedIndices = zeros(length(allVocalizedPairs), 1);
            
            %%- 01: BUILD WORD PAIRS
            % get word pairs in this session-block
            targetWords = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
            targetWords = {targetWords(3:end).name};
            wordPairs = createWithinVocalizedWordGroups(targetWords)
            
            %%- BUILD FEATURE MATRIX FOR EACH WORD PAIR
            sessionBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});

            %%- set plotting directories and meta data output
            figureDir = strcat('./Figures/', subj, '/reinstatement/within_blocks_vocalizationWord/');
            matDir = strcat('./Figures/', subj, '/reinstatement_mat/within_blocks_vocalizationWord/');
            figureFile = strcat(figureDir, sessions{iSesh}, '-', num2str(blocks{iBlock}));
            matFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}));  
            if ~exist(figureDir)
                mkdir(figureDir)
            end
            if ~exist(matDir)
                mkdir(matDir)
            end
            
            wordpairs = [wordPairs{:}];
            logFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), '.txt');
            fid = fopen(logFile, 'w');
            fprintf(fid, '%6s \n', 'Block(i) Target Words (vocalized):');
            fprintf(fid, '%6s \n', targetWords{:}); 
            fprintf(fid, '%6s \n', 'Word Pairs:');
            fprintf(fid, '%6s \n', wordpairs{:}); 
            
            fig = figure;
            clim = [0 0]; %initialize colorbar
            
            fa = {};
            %%- loop through every unique word pair of vocalized words
            for iWord=1:length(wordPairs),
                wordone = wordPairs{iWord}{1}; % first vocalized word
                wordtwo = wordPairs{iWord}{2}; % second vocalized word

                %%- 02: BUILD FEATURE MATRICES
                %%- get pair feature matrices for every vocalized word
                %%- comparison
                [pairFeatureMat1, pairFeatureMat2] = buildWithinPairFeatureMat(wordone, wordtwo, sessionBlockDir);
                
                size(pairFeatureMat1)
                size(pairFeatureMat2)
                
                % change time/features dimensions
                pairFeatureMat1 = permute(pairFeatureMat1, [1 3 2]);
                pairFeatureMat2 = permute(pairFeatureMat2, [1 3 2]);
                
                %%- 03: BUILD REINSTATEMENT MATRICES
                [eventRein, featureRein] = compute_reinstatement(pairFeatureMat1, pairFeatureMat2);
                size(eventRein)
                size(featureRein)
                
                timeZero = 16;
                ticks = [6:10:56]; % for 6 seconds of data
                labels = [-1:1:4];
                LT = 1.5; % line thickness
                
                %%- 04: PLOTTING
                fig;
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
                
                % save the relevant mat files
                save(strcat(matFile, '|', wordPairs{iWord}{1}, '_', wordPairs{iWord}{2}, '.mat'), 'eventRein', 'featureRein');
                
                checkOne = strjoin({wordone, wordtwo}, '_');
                checkTwo = strjoin({wordtwo, wordone}, '_');
                checkOne
                checkTwo
                
                %%- Check if this pair was already analyzed
                if (ismember(checkOne, allVocalizedPairs) ||...
                    ismember(checkTwo, allVocalizedPairs))
%                     
%                     index = cellfun(@(x) strcmp(checkOne, x), allVocalizedPairs, 'UniformOutput', 0);
%                     if isempty(find([index{:}] == 1))
%                         index = cellfun(@(x) strcmp(checkOne, x), allVocalizedPairs, 'UniformOutput', 0);
%                     end
%                     index = find([index{:}] == 1);
%                     
% %                     allVocalizedPairs
% %                     index
%                     allVocalizedPairs{index};
%                     allVocalizedIndices(index) = 1;
                end
            end % end of loop through words            
            
            % change the color limit to the max in the group for comparison
            for i=1:length(fa)
                fa{i}.CLim = clim;
            end
            
            % change figure dimensions before saving
            fig = gcf;
            fig.PaperUnits = 'inches';
            pos = [5.1667 0.6806 9.9722 10.3889];
            fig.PaperPosition = pos;
            
            %%- Save the image
            print(figureFile, '-dpng', '-r0')
            savefig(figureFile)
        end
    end
% end