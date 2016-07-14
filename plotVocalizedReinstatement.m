close all
clc
clear all


%%- Section for Pulling reinstatement matrices already produced and just
%%make different plots
ANALYSIS_TYPE = {'within_blocks', 'across_blocks'};
EVENT_SYNC = {'probeon', 'vocalizationWord'};

subj = 'NIH034';
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
% if strcmp(subj, 'NIH039')
%     sessions = sessions([1,2,4]);
% elseif strcmp(subj, 'NIH034')
%     sessions = sessions([3, 4]);
% end

%%- Parameter settings
subj = 'NIH034';
timeLock = 'vocalization';
referenceType = 'bipolar';
typeTransform = 'morlet';
typeAnalysis = 'across_blocks_vocalizationWord';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataDir = strcat('./Figures/', subj);
dataDir = fullfile(dataDir, 'reinstatement_mat', TYPE_TRANSFORM, typeAnalysis);
matFiles = dir(strcat(dataDir, '/*.mat')); % get all the mat files for this analysis
matFiles = {matFiles(:).name};

% initialize
allVocalizedPairs = {'CLOCK_JUICE', 'CLOCK_PANTS', 'CLOCK_BRICK', 'CLOCK_GLASS', 'CLOCK_CLOCK',...
                    'BRICK_JUICE', 'BRICK_PANTS', 'BRICK_BRICK', 'BRICK_GLASS', ...
                    'PANTS_JUICE', 'PANTS_PANTS', 'PANTS_GLASS', 'GLASS_JUICE', ...
                    'GLASS_GLASS', 'JUICE_JUICE'};
allVocalizedIndices = zeros(15, 1);
avgeEventReinMat = cell(15,1);
avgeFeatureReinMat = cell(15,1);

%%- average results across the entire subject tests
for iFile=1:length(matFiles)
    data = load(fullfile(dataDir, matFiles{iFile}));
    
    %%- extract event and feature reinstatements
    tempVocalizedIndices = data.allVocalizedIndices;
    eventReinMat = data.eventReinMat;
    featureReinMat = data.featureReinMat;
    
    buffIndex = 1;
    for iVocal=1:length(tempVocalizedIndices)
        if tempVocalizedIndices(iVocal) == 1 % this target pair is stored in this file
            if isempty(avgeEventReinMat{iVocal})
                avgeEventReinMat{iVocal} = eventReinMat{buffIndex};
                avgeFeatureReinMat{iVocal} = featureReinMat{buffIndex};
            else   
                avgeEventReinMat{iVocal} = cat(1, avgeEventReinMat{iVocal}, eventReinMat{buffIndex});
                avgeFeatureReinMat{iVocal} = cat(1, avgeFeatureReinMat{iVocal}, featureReinMat{buffIndex});
            end
            
            % increment index for reinstatement maps
            buffIndex = buffIndex + 1;
        end
    end
end

size(avgeEventReinMat)

%%- 04: PLOTTING
fig = figure;
clim = [0 0]; %initialize colorbar
fa = {};
LT = 1.5 %line thickness

if strcmp(typeTransform, 'morlet')
    ticks = [6:10:46];
    labels = [-1:1:3];
    timeZero = 16;
else
    timeZero = 5;
    ticks = [1:2:11];
    labels = [-2:1:3];
end


%%- Plot each word pairing separately
for iPlot=1:length(allVocalizedPairs)
    eventRein = avgeEventReinMat{iPlot};
    featureRein = avgeFeatureReinMat{iPlot};

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
fig.PaperUnits = 'inches';
pos = [0    0.6806   20.0000   10.2917];
fig.PaperPosition = pos;

%%- Save the image
dataDir = strcat('./Figures/', subj);
dataDir = fullfile(dataDir, 'reinstatement', TYPE_TRANSFORM, typeAnalysis);
figureFile = fullfile(dataDir, 'subject_summary');

print(figureFile, '-dpng', '-r0')
savefig(figureFile)
