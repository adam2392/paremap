%         -- Creating matrices that are eventsXtimeXfeatures, which can be
%         fed into compute_reinstatement.m
%         -- This is done for ACROSS blocks analysis of the paremap task

clear all;
close all;
clc;

%% PARAMETERS FOR RUNNING PREPROCESS
subj = 'NIH034';
sessNum = [0, 1, 2];

addpath('./m_reinstatement/');
%% LOAD EVENTS STRUCT AND SET DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load events and set behavioral directories                   ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
% eegRootDirHome = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
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
TYPE_TRANSF = 'morlet_spec';
dataDir = strcat('condensed_data_', subj);
dataDir = fullfile(dataDir, TYPE_TRANSF);
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
sessions = sessions(3:end);
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

%%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:length(blocks)-1,
        % get word pairs in this session-block (i)
        firstwordpairs = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        firstwordpairs = {firstwordpairs(3:end).name};
        % get the word pairs in the next block (i+1)
        secondwordpairs = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock+1}));
        secondwordpairs = {secondwordpairs(3:end).name};
        
        firstwordpairs
        secondwordpairs
        
        %%- CREATE ACROSS WORD PAIRS GROUPS 
        [sameWordGroup, reverseWordGroup, probeWordGroup, targetWordGroup, diffWordGroup] = createAcrossWordGroups(firstwordpairs, secondwordpairs);
        
        % sessionblock directories for block(i) and block(i+1)
        sessionFirstBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        sessionSecondBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock+1});
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
        
        %%- Build Similarity Matrics
        % same Pairs
        [eventSame, featureSame] = compute_reinstatement(samePairFeatureMat1, samePairFeatureMat2);
        [eventDiff, featureDiff] = compute_reinstatement(diffPairFeatureMat1, diffPairFeatureMat2);
        [eventReverse, featureReverse] = compute_reinstatement(reversePairFeatureMat1, reversePairFeatureMat2);
        [eventProbe, featureProbe] = compute_reinstatement(probePairFeatureMat1, probePairFeatureMat2);
        [eventTarget, featureTarget] = compute_reinstatement(targetPairFeatureMat1, targetPairFeatureMat2);
        
        size(squeeze(mean(eventDiff(:, :, :),1)))
        
        % rand sample down the different word pair feature mat -> match
        % size
        randIndices = randsample(size(eventDiff,1), size(eventSame,1));
        tempDiff = eventDiff(randIndices,:,:);
        
        % set linethickness
        LT = 1.5;
        
        %%- Plotting
        figure
        subplot(321)
        imagesc(squeeze(mean(eventSame(:, :, :),1)));
        title(['Same Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(tempDiff,1)), ' events'])
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = [0:10:55];
        ax.YTickLabel = [-1:1:5];
        ax.XTick = [0:10:55];
        ax.XTickLabel = [-1:1:4];
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        colorbar();
        clim = get(gca, 'clim');
        hold on
        plot(get(gca, 'xlim'), [10 10], 'k', 'LineWidth', LT)
        plot([10 10], get(gca, 'ylim'), 'k', 'LineWidth', LT)
       
        subplot(323);
        imagesc(squeeze(mean(tempDiff(:, :, :),1)));
        title(['Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1), ...
            ' vs ', num2str(iBlock)])
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = [0:10:55];
        ax.YTickLabel = [-1:1:5];
        ax.XTick = [0:10:55];
        ax.XTickLabel = [-1:1:4];
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        colorbar();
        set(gca, 'clim', clim);
        hold on
        plot(get(gca, 'xlim'), [10 10], 'k', 'LineWidth', LT)
        plot([10 10], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        subplot(325);
        imagesc(squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(tempDiff(:, :, :),1)));
        title(['Same-Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1), ...
            ' vs ', num2str(iBlock)])
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = [0:10:55];
        ax.YTickLabel = [-1:1:5];
        ax.XTick = [0:10:55];
        ax.XTickLabel = [-1:1:4];
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        colorbar();
        hold on
        plot(get(gca, 'xlim'), [10 10], 'k', 'LineWidth', LT)
        plot([10 10], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        %%- reverse, probe, target
        subplot(322)
        imagesc(squeeze(mean(eventReverse(:, :, :),1)));
        title(['Reverse Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(eventReverse,1)), ' events'])
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = [0:10:55];
        ax.YTickLabel = [-1:1:5];
        ax.XTick = [0:10:55];
        ax.XTickLabel = [-1:1:4];
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        colorbar();
        set(gca, 'clim', clim);
        hold on
        plot(get(gca, 'xlim'), [10 10], 'k', 'LineWidth', LT)
        plot([10 10], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        subplot(324)
        imagesc(squeeze(mean(eventProbe(:, :, :),1)));
        title(['Probe Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(eventProbe,1)), ' events'])
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = [0:10:55];
        ax.YTickLabel = [-1:1:5];
        ax.XTick = [0:10:55];
        ax.XTickLabel = [-1:1:4];
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        colorbar();
        set(gca, 'clim', clim);
        hold on
        plot(get(gca, 'xlim'), [10 10], 'k', 'LineWidth', LT)
        plot([10 10], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        subplot(326)
        imagesc(squeeze(mean(eventTarget(:, :, :),1)));
        title(['Target Pairs Cosine Similarity for Block ', num2str(iBlock-1), ' vs ',...
            num2str(iBlock), ' with ', num2str(size(eventTarget,1)), ' events'])
        hold on
        xlabel('Time (seconds)');
        ylabel('Time (seconds)');
        ax = gca;
        axis square
        ax.YTick = [0:10:55];
        ax.YTickLabel = [-1:1:5];
        ax.XTick = [0:10:55];
        ax.XTickLabel = [-1:1:4];
        colormap('jet');
        set(gca,'tickdir','out','YDir','normal');
        set(gca, 'box', 'off');
        colorbar();
        set(gca, 'clim', clim);
        hold on
        plot(get(gca, 'xlim'), [10 10], 'k', 'LineWidth', LT)
        plot([10 10], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        %%- Save Image
        figureDir = strcat('./Figures/', subj, '/reinstatement/across_blocks/');
        figureFile = strcat(figureDir, sessions{iSesh}, '-', num2str(blocks{iBlock}), 'vs',num2str(blocks{iBlock+1}));
        if ~exist(figureDir)
            mkdir(figureDir)
        end
        saveas(gca, figureFile, 'png')
        savefig(figureFile)
        
        %%- save reinstatement matrices
        save(strcat(figureFile, '.mat'), 'eventSame', 'featureSame', ...
                                        'eventReverse', 'featureReverse', ...
                                        'eventDiff', 'featureDiff', ...
                                        'eventTarget', 'featureTarget', ...
                                        'eventProbe', 'featureProbe');
        
        pause(0.1);
    end
end