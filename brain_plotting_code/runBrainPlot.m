clear all;
clc;
close all;

disp('runnning')
% brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'delta');
% brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'theta');
% brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'alpha');
% brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'beta');
% % % brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'low gamma');
% brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'high gamma');
brainPlotSpectralPower('NIH034', 'morlet', 'bipolar', 'HFO');

disp('nih037')
% brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'delta');
% brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'theta');
% brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'alpha');
% brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'beta');
% brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'low gamma');
% brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'high gamma');
% brainPlotSpectralPower('NIH037', 'morlet', 'bipolar', 'HFO');

% brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'delta');
% brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'theta');
% brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'alpha');
% brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'beta');
% brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'low gamma');
% brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'high gamma');
brainPlotSpectralPower('NIH039', 'morlet', 'bipolar', 'HFO');

% save video
subj = 'NIH034';
typeTransform = 'morlet';
referenceType = 'bipolar';
frequencyBand = 'HFO';

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

% figure Dirs
brainFigDir = fullfile(eegRootDir, '/brain_plotting_code/Figures/', ...
    strcat(typeTransform, '_', referenceType, '_', frequencyBand, '_within_blocks'));

groupOneDir = fullfile(brainFigDir, 'groupone');
brainPlots = dir(fullfile(groupOneDir, '*.png'));
brainPlots = {brainPlots.name};

groupOneOutput = fullfile(brainFigDir, strcat(subj,'-', frequencyBand, '-groupone.avi'));
delete(groupOneOutput);
outputVideo = VideoWriter(groupOneOutput);
outputVideo.FrameRate = 1;
outputVideo.Quality = 100;
open(outputVideo);
for iPlot=1:length(brainPlots)
    figFile = fullfile(groupOneDir, brainPlots{iPlot});
    
    img = imread(figFile);
    disp(['Writing image ', iPlot, ' to video...'])
    writeVideo(outputVideo, img);
end
close(outputVideo);
disp('Done!');

groupTwoDir = fullfile(brainFigDir, 'grouptwo');
brainPlots = dir(fullfile(groupTwoDir, '*.png'));
brainPlots = {brainPlots.name};

groupOneOutput = fullfile(brainFigDir, strcat(subj,'-', frequencyBand, '-grouptwo.avi'));
delete(groupOneOutput);
outputVideo = VideoWriter(groupOneOutput);
outputVideo.FrameRate = 1;
outputVideo.Quality = 100;
open(outputVideo);
for iPlot=1:length(brainPlots)
    figFile = fullfile(groupTwoDir, brainPlots{iPlot});
    
    img = imread(figFile);
    disp(['Writing image ', iPlot, ' to video...'])
    writeVideo(outputVideo, img);
end
close(outputVideo);
disp('Done!');