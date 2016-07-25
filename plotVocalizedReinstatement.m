function plotVocalizedReinstatement(subj, typeTransform, referenceType, typeAnalysis)
close all
clc

%%- Parameter settings
% subj = 'NIH034';
% referenceType = 'bipolar';
% typeTransform = 'morlet';
% typeAnalysis = 'within_blocks_vocalizationWord';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
eegRootDirHome = '/Volumes/NIL_PASS';
eegRootDirJhu = '/home/adamli/paremap';
% Determine which directory we're working with automatically
if     ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

dataDir = strcat('./Figures/', subj);
analysisDir = strcat(typeTransform, '_', referenceType, '_', typeAnalysis);
dataDir = fullfile(eegRootDir, dataDir, 'reinstatement_mat', analysisDir);
matFiles = dir(strcat(dataDir, '/*.mat')); % get all the mat files for this analysis
matFiles = {matFiles(:).name}

% initialize
allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                    'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                    'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', 'GLASS_JUICE', ...
                    'GLASS_GLASS', 'JUICE_JUICE'};
allVocalizedIndices = zeros(15, 1);
avgeEventReinMat = cell(15,1);
% avgeFeatureReinMat = cell(15,1);

%%- average results across the entire subject tests
for iFile=1:length(matFiles)
    data = load(fullfile(dataDir, matFiles{iFile}));
    
    %%- extract event and feature reinstatements
    tempVocalizedIndices = data.allVocalizedIndices;
    eventReinMat = data.eventReinMat;
%     featureReinMat = data.featureReinMat;
    
    % loop through each index of targetword comparisons
    for iVocal=1:length(tempVocalizedIndices)
        if tempVocalizedIndices(iVocal) == 1 % this target pair is stored in this file
            if isempty(avgeEventReinMat{iVocal})
                avgeEventReinMat{iVocal} = eventReinMat{iVocal};
%                 avgeFeatureReinMat{iVocal} = featureReinMat{iVocal};
            else   
                avgeEventReinMat{iVocal} = cat(1, avgeEventReinMat{iVocal}, eventReinMat{iVocal});
%                 avgeFeatureReinMat{iVocal} = cat(1, avgeFeatureReinMat{iVocal}, featureReinMat{iVocal});
            end
        end
    end
end

% get the min number of events, so everything is same # of events
size(avgeEventReinMat)
for iCell=1:length(avgeEventReinMat)
    if iCell==1
        minNumEvents = size(avgeEventReinMat{iCell}, 1);
    else
        minNumEvents = min(minNumEvents, size(avgeEventReinMat{iCell}, 1));
    end
end
% downsample all events from targetword comparisons
for iCell=1:length(avgeEventReinMat)
    randIndices = randsample(size(avgeEventReinMat{iCell},1), minNumEvents);
    events = avgeEventReinMat{iCell};
    events = events(randIndices,:,:);
    
    avgeEventReinMat{iCell} = events;
end

%%- 04: PLOTTING
fig = figure;
clim = [0 0]; %initialize colorbar
fa = {};
LT = 1.5 %line thickness

% load in an example file to get the -> labels, ticks and timeZero
dataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(dataDir, strcat(typeTransform, '_', referenceType,  '_vocalization_sessiontargetwords'));
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

pairDirs = dir(fullfile(dataDir, sessions{1}, blocks{1}));
exampleDir = fullfile(dataDir, sessions{1}, blocks{1}, pairDirs(4).name);
channelData = dir(exampleDir);
data = load(fullfile(exampleDir, channelData(4).name));
data = data.data;
timeTicks = data.waveT(:,2);

ticks = 1:5:length(timeTicks);
labels = timeTicks(ticks);
timeZero = data.timeZero;

%%- Plot each word pairing separately
for iPlot=1:length(allVocalizedPairs)
    eventRein = avgeEventReinMat{iPlot};
%     featureRein = avgeFeatureReinMat{iPlot};

    % split words
    wordSplit = strsplit(allVocalizedPairs{iPlot}, '_');
    wordone = wordSplit{1};
    wordtwo = wordSplit{2};

    fa{iPlot} = subplot(4, 4, iPlot);
    imagesc(squeeze(mean(eventRein(:,:,:),1)));
    title({['Cosine Similarity for Subject ', subj, ' for '], [wordone, ' vs ', wordtwo, ...
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
end
% change the color limit to the max in the group for comparison
for i=1:length(allVocalizedPairs)
    if allVocalizedIndices(i) == 1
        fa{i}.CLim = clim;
    end
end

% change figure dimensions before saving
fig = gcf;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
pos = [0    0.6667   17.5972   10.4028];
fig.PaperPosition = pos;

%%- Save the image
dataDir = strcat('./Figures/', subj);
dataDir = fullfile(dataDir, 'reinstatement', analysisDir);
figureFile = fullfile(dataDir, 'subject_summary');

print(figureFile, '-dpng', '-r0')
savefig(figureFile)
end
