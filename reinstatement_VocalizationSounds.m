function reinstatement_VocalizationSounds(subj, typeTransform, referenceType)
% subj = 'NIH034';
% typeTransform = 'morlet';
% referenceType = 'bipolar';
close all;

addpath('./reinstatement_vocalizedWords/');
% determine which directory path to use
eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
eegRootDirVolume = '/Volumes/NIL_PASS';
eegRootDirJhu = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
elseif ~isempty(dir(eegRootDirVolume)), eegRootDir = eegRootDirVolume;
elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjDataDir = strcat('/condensed_data_', subj);
dataDir = fullfile(eegRootDir, subjDataDir, strcat(typeTransform, '_', referenceType, '_', 'vocalization_sessiontargetwords'));
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

figureDir = fullfile(eegRootDir, '/Figures/', subj, 'reinstatement_sounds',...
    strcat(typeTransform, '_', referenceType, '_within_blocks_vocalizationWord/'));
matDir = fullfile(eegRootDir, '/Figures/', subj, 'reinstatement_mat_sounds',...
    strcat(typeTransform, '_', referenceType, '_within_blocks_vocalizationWord/'));
if ~exist(figureDir, 'dir') mkdir(figureDir); end
if ~exist(matDir, 'dir')    mkdir(matDir);    end
% load in an example file to get the labels, ticks and timeZero
pairDirs = dir(fullfile(dataDir, sessions{1}, blocks{1}));
exampleDir = fullfile(dataDir, sessions{1}, blocks{1}, pairDirs(4).name);
channelData = dir(exampleDir);
data = load(fullfile(exampleDir, channelData(4).name));
data = data.data;
timeTicks = data.waveT(:,2);

ticks = 6:10:length(timeTicks);
labels = timeTicks(ticks);
timeZero = data.timeZero;

disp(['Analyzing ', sessions{:}, ' for ', subj]);

group_one = {'CLOCK', 'BRICK'};
group_two = {'PANTS', 'GLASS', 'JUICE'};
LT = 1.5; % line thickness

%% RUN ANALYSIS PER SESSION-BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Separate into two groups and compare reinstatement        -----%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSesh=1:length(sessions),
    for iBlock=1:length(blocks),
        fprintf('%6s \n', strcat('On session ', num2str(iSesh), ' and block ', num2str(iBlock)));
        
        % get word pairs in this session-block
        targetWords = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        targetWords = {targetWords(3:end).name};
        sessionBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        
        targetOne = {};
        targetTwo = {};
        
        % concatenated Z-scored power matrices
        featureMatGroupOne = [];
        featureMatGroupTwo = [];
        %%- LOOP THROUGH TARGET WORDS AND PUT INTO 1 OF 2 GROUPS
        for iWord=1:length(targetWords)
            featureMat = buildFeatureMat(targetWords{iWord}, sessionBlockDir);
            if find(strcmp(group_one, targetWords{iWord}))
                targetOne{end+1} = targetWords{iWord};
                
                if isempty(featureMatGroupOne)
                    featureMatGroupOne = featureMat;
                else
                    featureMatGroupOne = cat(1, featureMatGroupOne, featureMat);
                end
            else
                targetTwo{end+1} = targetWords{iWord};
                
                if isempty(featureMatGroupTwo)
                    featureMatGroupTwo = featureMat;
                else
                    featureMatGroupTwo = cat(1, featureMatGroupTwo, featureMat);
                end
            end
        end
        
        %%- CREATE A PAIRWISE COMPARISON OF ALL PAIRS NOT IN SAME GROUPS
        [same1reinstatement_inputone, same1reinstatement_inputtwo] = buildReinstatementInput(featureMatGroupOne, featureMatGroupOne, 1);
        [same2reinstatement_inputone, same2reinstatement_inputtwo] = buildReinstatementInput(featureMatGroupTwo, featureMatGroupTwo, 1);
        [diffreinstatement_inputone, diffreinstatement_inputtwo] = buildReinstatementInput(featureMatGroupOne, featureMatGroupTwo, 0);
        
        randIndices = randsample(size(diffreinstatement_inputone, 1), size(same1reinstatement_inputone, 1));
        diffreinstatement_inputone = diffreinstatement_inputone(randIndices, :, :);
        diffreinstatement_inputtwo = diffreinstatement_inputtwo(randIndices, :, :);

        [eventSame1, featureSame1] = compute_reinstatement(same1reinstatement_inputone, same1reinstatement_inputtwo);
        [eventSame2, featureSame2] = compute_reinstatement(same2reinstatement_inputone, same2reinstatement_inputtwo);
        [eventDiff, featureDiff] = compute_reinstatement(diffreinstatement_inputone, diffreinstatement_inputtwo);
        
        %%- save01: mat files
        matFile = fullfile(matDir, strcat(sessions{iSesh}, '-', blocks{iBlock}, '.mat'));
        save(matFile, 'eventSame1', 'featureSame1', ...
            'eventSame2', 'featureSame2', ...
            'eventDiff', 'featureDiff',...
            'group_one', 'group_two');
        
        %%- save02: reinstatement plots
        clim = [3 -4];
        figure;
        FA = {};
        FA{1} = subplot(311);
        imagesc(squeeze(mean(eventSame1(:, :, :),1)));
        title(['Same Sound Comparisons of CLOCK/BRICK for ', sessions{iSesh},...
            ' and ', num2str(iBlock-1)])
        colorbar(); colormap('jet'); 
        tempclim = get(gca, 'clim');
        clim(1) = min(tempclim(1), clim(1));
        clim(2) = max(tempclim(2), clim(2));
        axis square; hold on
        xlabel('Time (seconds)'); ylabel('Time (seconds)');
        ax = gca;
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        set(gca,'tickdir','out','YDir','normal', 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT);
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT);
        
        FA{2} = subplot(312);
        imagesc(squeeze(mean(eventSame2(:, :, :),1)));
        title(['Same Sound Comparisons of GLASS/JUICE/PANTS for ', sessions{iSesh},...
            ' and ', num2str(iBlock-1)])
        colorbar(); colormap('jet'); 
        tempclim = get(gca, 'clim');
        clim(1) = min(tempclim(1), clim(1));
        clim(2) = max(tempclim(2), clim(2));
        axis square; hold on
        xlabel('Time (seconds)'); ylabel('Time (seconds)');
        ax = gca;
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        set(gca,'tickdir','out','YDir','normal', 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT);
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT);
        
        FA{3} = subplot(313);
        imagesc(squeeze(mean(eventDiff(:, :, :),1)));
        title(['Diff Sound Comparisons for ', sessions{iSesh},...
            ' and ', num2str(iBlock-1)])
        colorbar(); colormap('jet'); 
        tempclim = get(gca, 'clim');
        clim(1) = min(tempclim(1), clim(1));
        clim(2) = max(tempclim(2), clim(2));
        axis square; hold on
        xlabel('Time (seconds)'); ylabel('Time (seconds)');
        ax = gca;
        ax.YTick = ticks;
        ax.YTickLabel = labels;
        ax.XTick = ticks;
        ax.XTickLabel = labels;
        set(gca,'tickdir','out','YDir','normal', 'box', 'off');
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT);
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT);
        
        % change the color limit to the max in the group for comparison
        for i=1:length(FA)
            FA{i}.CLim = clim;
        end
        
        % change figure dimensions before saving
        fig = gcf;
        fig.Units = 'inches';
        fig.PaperUnits = 'inches';
        pos = [0    0.6667   17.5972   10.4028];
        fig.PaperPosition = pos;

        %%- Save the image
        figureFile = fullfile(figureDir, strcat(sessions{iSesh}, '-', blocks{iBlock}));
        print(figureFile, '-dpng', '-r0')
        
        clear eventSame1 eventSame2 eventDiff featureDiff featureSame1 featureSame2
    end
end



end

