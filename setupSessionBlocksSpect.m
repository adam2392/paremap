function setupSessionBlocksSpect(subj, typeTransform, referenceType, blocksComp)
clc;
close all;

% subj = 'NIH034';
% typeTransform = 'multitaper';
% referenceType = 'bipolar';
% blocksComp = 'within_blocks';

%% PARAMETERS FOR RUNNING PREPROCESS
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
%% LOAD PREPROCESSED DATA DIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ STEP 2: Load data from Dir and create eventsXfeaturesxTime    ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjDataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(subjDataDir, strcat(typeTransform, '_', referenceType, '_', 'vocalization_sessiontargetwords'));
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};
sessions

% all target word comparisons that exist 
allVocalizedWords = {'CLOCK', 'JUICE', 'PANTS', 'BRICK', 'GLASS'};

%%- SAVING FIGURES OPTIONS
figureDir = strcat('./Figures/', subj, '/reinstatement/', typeTransform, '_', referenceType, '_', 'within_blocks_vocalizationWord/');
matDir = strcat('./Figures/', subj, '/reinstatement_mat/', typeTransform, '_', referenceType, '_', 'within_blocks_vocalizationWord/');
if ~exist(figureDir, 'dir') mkdir(figureDir); end
if ~exist(matDir, 'dir')    mkdir(matDir);    end

LT = 1.5; % line thickness

% load in an example file to get the labels, ticks and timeZero
pairDirs = dir(fullfile(dataDir, sessions{1}, blocks{1}));
exampleDir = fullfile(dataDir, sessions{1}, blocks{1}, pairDirs(4).name);
channelData = dir(exampleDir);
data = load(fullfile(exampleDir, channelData(4).name));
data = data.data;
timeTicks = data.waveT(:,2);

ticks = 1:5:length(timeTicks);
labels = timeTicks(ticks);
timeZero = data.timeZero;
  
%% CREATE VOCALIZED WORD GROUPS
%%- LOOP THROUGH SESSIONS AND BLOCKS
%%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    sessionPowerMat = cell(5, 1); % per session power matrix
    
    if ~exist(strcat(matDir, '/channels_vocalized_session/'))
        mkdir(strcat(matDir, '/channels_vocalized_session/'));
    end
    matSessionFile = strcat(matDir, '/channels_vocalized_session/', 'spect_', sessions{iSesh});
    
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:length(blocks),
        fprintf('%6s \n', strcat('On session ', num2str(iSesh), ' and block ', num2str(iBlock)));
        
        % get word pairs in this session-block
        targetWords = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        targetWords = {targetWords(3:end).name};
        
        %%- set plotting directories and meta data output
        matBlockFile = strcat(matDir, '/channels_vocalized_blocks/', 'spect_', sessions{iSesh}, '_', blocks{iBlock});
        if ~exist(strcat(matDir, '/channels_vocalized_blocks/'))
            mkdir(strcat(matDir, '/channels_vocalized_blocks/'));
        end
        
        sessionBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        blockPowerMat = cell(5,1); 
        
        %%- Load in power Data and average
        for iTarget=1:length(targetWords)
            currentTarget = targetWords{iTarget};
            wordDir = fullfile(sessionBlockDir, currentTarget);
            chanFiles = dir(wordDir);
            chanFiles = {chanFiles(3:end).name};
            
            avgePowerMat = [];
              
            for iChan=1:length(chanFiles)
                filePath = fullfile(wordDir, chanFiles{iChan});
                data = load(filePath);
                data = data.data;
                
                %%- Concatenate all frequency vectors into feature vector
                if isempty(avgePowerMat)
                    avgePowerMat = data.powerMatZ;
                else
                    avgePowerMat = cat(2, avgePowerMat, data.powerMatZ);
                end 
            end
            vocalizedIndex = find(strcmp(allVocalizedWords, currentTarget));
            blockPowerMat{vocalizedIndex} = avgePowerMat;
            
            if isempty(sessionPowerMat{vocalizedIndex})
                sessionPowerMat{vocalizedIndex} = avgePowerMat;
            else
                sessionPowerMat{vocalizedIndex} = cat(1, sessionPowerMat{vocalizedIndex}, avgePowerMat);
            end
        end
        
        %%- save per session-block and clear struct
        save(strcat(matBlockFile, '.mat'), 'blockPowerMat', 'allVocalizedWords');
    end
    
    %%- save mat per session
    save(strcat(matSessionFile, '.mat'), 'sessionPowerMat', 'allVocalizedWords');
end
end