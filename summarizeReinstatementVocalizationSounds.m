function summarizeReinstatementVocalizationSounds(subj, typeTransform, referenceType)
% subj = 'NIH039';
% typeTransform = 'multitaper';
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

LT = 1.5;
%%- LOOP THROUGH EACH MAT FILE AND AVERAGE TOGETHER
matFiles = dir(fullfile(matDir, '*.mat'));
matFiles = {matFiles.name};
beginindex = 1;
endindex = 0;
groupOneMat = zeros(1000*length(matFiles), length(timeTicks), length(timeTicks));
groupTwoMat = zeros(1000*length(matFiles), length(timeTicks), length(timeTicks));
groupDiffMat = zeros(1000*length(matFiles), length(timeTicks), length(timeTicks));

for iMat=1:length(matFiles)
    data = load(fullfile(matDir, matFiles{iMat}));
    eventSame1 = data.eventSame1;
    eventSame2 = data.eventSame2;
    eventDiff = data.eventDiff;
    
    if size(eventSame1,1) ~= size(eventSame2,1) || size(eventSame2,1) ~= size(eventDiff,1)
        disp('ERROR');
        pause;
    end
    
    endindex = endindex + size(eventSame1,1);
    groupOneMat(beginindex:endindex,:,:) = eventSame1;
    groupTwoMat(beginindex:endindex,:,:) = eventSame2;
    groupDiffMat(beginindex:endindex,:,:) = eventDiff;
    beginindex = endindex+1;
end

groupOneMat(beginindex:end,:,:) = [];
groupTwoMat(beginindex:end,:,:) = [];
groupDiffMat(beginindex:end,:,:) = [];

% downsample all to same # of events
minsamples = min([size(groupOneMat,1), size(groupTwoMat,1), size(groupDiffMat,1)]);
randindices = randsample(size(groupOneMat,1), minsamples);
groupOneMat = groupOneMat(randindices,:,:);
randindices = randsample(size(groupTwoMat,1), minsamples);
groupTwoMat = groupTwoMat(randindices,:,:);
randindices = randsample(size(groupDiffMat,1), minsamples);
groupDiffMat = groupDiffMat(randindices,:,:);

clim = [3 -4];
FA = {};
figure;
FA{1} = subplot(311);
imagesc(squeeze(mean(groupOneMat,1)));
title(['BRICK/CLOCK Sounds ', subj, ' summary'])
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
imagesc(squeeze(mean(groupTwoMat,1)));
title(['GLASS/JUICE/PANTS Sounds ', subj, ' summary'])
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
imagesc(squeeze(mean(groupTwoMat,1) - mean(groupOneMat,1)));
title(['Differences Between Reinstatement Sounds ', subj, ' summary'])
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
for i=1:length(FA)-1
    FA{i}.CLim = clim;
end

% change figure dimensions before saving
fig = gcf;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
pos = [0    0.6667   17.5972   10.4028];
fig.PaperPosition = pos;

%%- Save the image
figureFile = fullfile(figureDir, strcat(subj, '-summary'));
print(figureFile, '-dpng', '-r0')
end
        