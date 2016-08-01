function brainPlotSpectralPower(subj, typeTransform, referenceType, frequencyBand)
subj = 'NIH037';
typeTransform = 'morlet';
referenceType = 'bipolar';
frequencyBand = 'HFO';
radius = 12.5;
subj
typeTransform
frequencyBand

% array of frequency bands
freqBandAr(1).name    = 'delta';
freqBandAr(1).rangeF  = [2 4];          %[2 4]
freqBandAr(2).name    = 'theta';
freqBandAr(2).rangeF  = [4 8];          %[4 8]
freqBandAr(3).name    = 'alpha';
freqBandAr(3).rangeF  = [8 16];         %[8 12]
freqBandAr(4).name    = 'beta';
freqBandAr(4).rangeF  = [16 32];        %[12 30]
freqBandAr(5).name    = 'low gamma';
freqBandAr(5).rangeF  = [32 80];        %[30 70]
freqBandAr(6).name    = 'high gamma';
freqBandAr(6).rangeF  = [80 160];       %[70 150]
freqBandAr(7).name    = 'HFO';
freqBandAr(7).rangeF  = [160 400];      %[150 400]
freqIndice = find(strcmp({freqBandAr.name}, frequencyBand));

load ROI.mat %load matrix of size number of ROIs X 3 containing the x,y,z, coordinates of each ROI

% determine which directory path to use
eegRootDirWork = '/Users/liaj/Documents/MATLAB/paremap';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home
eegRootDirVolume = '/Volumes/NIL_PASS';
eegRootDirJhu = '/home/adamli/paremap';

% Determine which directory we're working with automatically
if ~isempty(dir(eegRootDirVolume)), eegRootDir = eegRootDirVolume;
elseif ~isempty(dir(eegRootDirWork)), eegRootDir = eegRootDirWork;
elseif ~isempty(dir(eegRootDirHome)), eegRootDir = eegRootDirHome;
elseif ~isempty(dir(eegRootDirJhu)), eegRootDir = eegRootDirJhu;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 1: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjDataDir = strcat('/condensed_data_', subj);
dataDir = fullfile(eegRootDir, subjDataDir, strcat(typeTransform, '_', referenceType, '_', 'targetWords'));

% load in an example file to get the labels, ticks and timeZero
pairDirs = dir(fullfile(dataDir, 'BRICK'));
exampleDir = fullfile(dataDir, 'BRICK', pairDirs(4).name);
data = load(fullfile(exampleDir));
data = data.data;
timeTicks = data.waveT(:,2);

ticks = 1:5:length(timeTicks);
labels = timeTicks(ticks);
timeZero = data.timeZero;

els = load_electrode_info(subj, 1);
elec_locs = [[els.x]' [els.y]' [els.z]'];%Create matrix of size electrodes X 3 containing the x,y,z, coordinates of each electrode
n_elecs = numel(els);
elecToROI = create_el_to_roi_matrix(ROI, elec_locs, radius);
% [plots, brains] = plot3brains_base;

%% GENERATE DATA TO PLOT
group_one = {'CLOCK', 'BRICK'};
group_two = {'PANTS', 'GLASS', 'JUICE'};

% figure Dirs
brainFigDir = fullfile(eegRootDir, '/brain_plotting_code/Figures/', ...
    strcat(typeTransform, '_', referenceType, '_', frequencyBand, '_within_blocks'));
brainMatDir = fullfile(eegRootDir, '/brain_plotting_code/Mats/', ...
    strcat(typeTransform, '_', referenceType, '_', frequencyBand, '_within_blocks'));
if ~exist(brainFigDir, 'dir') mkdir(brainFigDir); end
if ~exist(brainMatDir, 'dir')    mkdir(brainMatDir);    end

fprintf('%6s \n', strcat('On subject ', subj));

% get word pairs in this session-block
targetWords = dir(fullfile(dataDir));
targetWords = {targetWords(3:end).name};
targetWordDir = fullfile(dataDir);

% concatenated Z-scored power matrices
featureMatGroupOne = [];
featureMatGroupTwo = [];
targetOne = {}; % clcok and brikc
targetTwo = {}; % pants, glass and juice
%%- LOOP THROUGH TARGET WORDS AND PUT INTO 1 OF 2 GROUPS
for iWord=1:length(targetWords)
    featureMat = buildBrainPlotFeatureMat(targetWords{iWord}, dataDir, freqIndice);
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

% average across events within this session-block
featureMatGroupOne = squeeze(mean(featureMatGroupOne, 1));
featureMatGroupTwo = squeeze(mean(featureMatGroupTwo, 1));

if strcmp(subj, 'NIH037')
    elecToROI = elecToROI(:, 1:end-1);
end

roi_vals_one = elecToROI * featureMatGroupOne; % roi values for each time point
roi_vals_two = elecToROI * featureMatGroupTwo; 

% save mat files
matFile = fullfile(brainMatDir, strcat(subj, '-', 'targetWordComp'));
save(matFile, 'roi_vals_one', 'roi_vals_two', 'targetOne', 'targetTwo');

%% CARRY OUT PLOTTING AND SAVE OF BRAIN FIGS
close all;
[plots1,brains1]=plot3brains_base(1);%This plots 3 brains and returns the handles to the axes and the brain surfaces
maxclim = max(max(roi_vals_one(:)), max(roi_vals_two(:)));
minclim = min(min(roi_vals_one(:)), min(roi_vals_two(:)));

tic;
h1 = [];
for iTime=12:size(featureMatGroupOne, 2)-9
    s1 = struct();
    s1.clim = [minclim, maxclim];

    use_rwb = 0;
    h1 = update3brains_v2(brains1,roi_vals_one(:, iTime),s1,...
        ['Looking at CLOCK/BRICK GROUP ', frequencyBand, ' for ', subj, ' at ', num2str(timeTicks(iTime))],'Z-scored Spectral Power',...
        use_rwb,h1);

    %%- save image
    fig = gcf;
    fig.Units = 'inches';
    fig.PaperUnits = 'inches';
    pos = [11.9375   -7.4896   19.3021   11.3750];
    fig.PaperPosition = pos;
    
    figureFile = fullfile(brainFigDir, 'groupone', strcat('groupone-',subj, '-', num2str(iTime)));
    print(figureFile, '-dpng', '-r0')
end
toc;

close all;
[plots2,brains2]=plot3brains_base(2);%This plots 3 brains and returns the handles to the axes and the brain surfaces  
tic;
h1 = [];
for iTime=12:size(featureMatGroupTwo, 2)-9
    s1 = struct();
    s1.clim = [minclim, maxclim];
    use_rwb = 0;
    h1 = update3brains_v2(brains2,roi_vals_two(:, iTime),s1,...
        ['Looking at GLASS/PANTS/JUICE GROUP ', frequencyBand, ...
        ' for ', subj, ' at ', num2str(timeTicks(iTime))],'Z-scored Spectral Power',...
        use_rwb,h1);

    %%- save image
    fig = gcf;
    fig.Units = 'inches';
    fig.PaperUnits = 'inches';
    pos = [11.9375   -7.4896   19.3021   11.3750];
    fig.PaperPosition = pos;
    
    figureFile = fullfile(brainFigDir, 'grouptwo',strcat('grouptwo-',subj, '-', num2str(iTime)));
    print(figureFile, '-dpng', '-r0')
end
toc;
end