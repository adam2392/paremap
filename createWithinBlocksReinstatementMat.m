%         -- Creating matrices that are eventsXtimeXfeatures, which can be
%         fed into compute_reinstatement.m
%         -- This is done for within blocks analysis of the paremap task
%        

function createWithinBlocksReinstatementMat(subj, VOCALIZATION, MATCHWORD)
close all;
% clear all;
% clc;

%% PARAMETERS FOR RUNNING PREPROCESS
if ~exist(subj)
    subj = 'NIH034';
end
sessNum = [0, 1, 2];
% if ~exist(VOCALIZATION)
%     VOCALIZATION = 0;
% end

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

%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if VOCALIZATION,
    TYPE_TRANSF = 'morlet_spec_vocalization';
elseif MATCHWORD
    TYPE_TRANSF = 'morlet_spec_matchword';
else
    TYPE_TRANSF = 'morlet_spec';
end

disp(TYPE_TRANSF)

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

%%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:length(blocks),
        % get word pairs in this session-block
        wordpairs = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        wordpairs = {wordpairs(3:end).name};
        
        %%- CREATE WITHIN WORD PAIRS GROUPS 
        [sameWordGroup, reverseWordGroup, diffWordGroup] = createWithinWordGroups(wordpairs);
        
        %%- BUILD FEATURE MATRIX 
        sessionBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        [samePairFeatureMat1, samePairFeatureMat2] = buildSamePairFeatureMat(sameWordGroup, sessionBlockDir);
        [reversePairFeatureMat1, reversePairFeatureMat2] = buildDiffPairFeatureMat(reverseWordGroup, sessionBlockDir);
        [diffPairFeatureMat1, diffPairFeatureMat2] = buildDiffPairFeatureMat(diffWordGroup, sessionBlockDir);

        size(samePairFeatureMat1)
        size(samePairFeatureMat2)
        size(diffPairFeatureMat1)
        size(diffPairFeatureMat2)
        
        % reshape matrix to events X time X features
        samePairFeatureMat1 = permute(samePairFeatureMat1, [1 3 2]);
        samePairFeatureMat2 = permute(samePairFeatureMat2, [1 3 2]);
        reversePairFeatureMat1 = permute(reversePairFeatureMat1, [1 3 2]);
        reversePairFeatureMat2 = permute(reversePairFeatureMat2, [1 3 2]);
        diffPairFeatureMat1 = permute(diffPairFeatureMat1, [1 3 2]);
        diffPairFeatureMat2 = permute(diffPairFeatureMat2, [1 3 2]);
        
        size(samePairFeatureMat1)
        size(samePairFeatureMat2)
        size(diffPairFeatureMat1)
        size(diffPairFeatureMat2)
        
        %%- Build Similarity Matrics
        % same Pairs
        [eventSame, featureSame] = compute_reinstatement(samePairFeatureMat1, samePairFeatureMat2);
        [eventReverse, featureReverse] = compute_reinstatement(reversePairFeatureMat1, reversePairFeatureMat2);
        [eventDiff, featureDiff] = compute_reinstatement(diffPairFeatureMat1, diffPairFeatureMat2);
        
        size(squeeze(mean(eventDiff(:, :, :),1)))
        
        % rand sample down the different word pair feature mat -> match
        % size
        randIndices = randsample(size(eventDiff,1), size(eventSame,1));
        eventDiff = eventDiff(randIndices,:,:);
        
        if VOCALIZATION,
            ticks = [6:10:56];
            labels = [-3:1:2];
            timeZero = 36;
        elseif MATCHWORD
            ticks = [6:10:56];
            labels = [-4:1:1];
            timeZero = 46;
        else
            ticks = [6:10:56];
            labels = [0:1:5];
            timeZero = 6;
        end
        
        % set linethickness
        LT = 1.5;
        
        %%- Save Image
        if VOCALIZATION,
            figureDir = strcat('./Figures/', subj, '/reinstatement/within_blocks_vocalization/');
            matDir = strcat('./Figures/', subj, '/reinstatement_mat/within_blocks_vocalization/');
        elseif MATCHWORD
            figureDir = strcat('./Figures/', subj, '/reinstatement/within_blocks_matchword/');
            matDir = strcat('./Figures/', subj, '/reinstatement_mat/within_blocks_matchword/');
        else
            figureDir = strcat('./Figures/', subj, '/reinstatement/within_blocks_probeon/');
            matDir = strcat('./Figures/', subj, '/reinstatement_mat/within_blocks_probeon/');
        end
        figureFile = strcat(figureDir, sessions{iSesh}, '-', num2str(blocks{iBlock}));
        matFile = strcat(matDir, sessions{iSesh}, '-', num2str(blocks{iBlock}));
        if ~exist(figureDir)
            mkdir(figureDir)
        end
        if ~exist(matDir)
            mkdir(matDir)
        end
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
        figure
        subplot(311)
        imagesc(squeeze(mean(eventSame(:, :, :),1)));
        title(['Same Pairs Cosine Similarity for Block ', num2str(iBlock-1) ...
            ' with ', num2str(size(eventDiff,1)), ' events'])
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
        clim = get(gca, 'clim');
        hold on
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
       
        subplot(312);
        imagesc(squeeze(mean(eventDiff(:, :, :),1)));
        title(['Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
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
        set(gca, 'clim', clim);
        hold on
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        subplot(313);
        imagesc(squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(eventDiff(:, :, :),1)));
        title(['Same-Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
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
        hold on
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        %%- SAve Imge
%         fig = gcf;
%         fig.PaperUnits = 'inches';
%         pos = [13.8750    4.5694   13.2639   14.0000];
%         fig.PaperPosition = pos;
        
        %%- Save Image
        print(figureFile, '-dpng', '-r0')
        savefig(figureFile)
        
        pause(0.1);
%         set(gca, 'clim', clim);  
    end % loop through blocks
end % loop through sessions
end
