clear all;
clc;

%% PARAMETERS FOR RUNNING PREPROCESS
subj = 'NIH034';
sessNum = [0, 1, 2];
DEBUG = 1;

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

sessions
blocks
events
dataDir

%%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:length(blocks),
        % get word pairs in this session-block
        wordpairs = dir(fullfile(dataDir, sessions{iSesh}, blocks{iBlock}));
        wordpairs = {wordpairs(3:end).name};
        %%- CREATE WITHIN WORD PAIRS GROUPS 
        [sameWordGroup, reverseWordGroup, diffWordGroup] = createWithinWordGroups(wordpairs);
        sessionBlockDir = fullfile(dataDir, sessions{iSesh}, blocks{iBlock});
        
        
        % loop through different words
        diffPairFeatureMat1 = [];
        diffPairFeatureMat2 = [];
        
        blockSpect = [];
        % loop through every word in the sameWordGroup
        for iWord=1:length(sameWordGroup),
            wordone = diffWordGroup{iWord}{1};
            wordtwo = diffWordGroup{iWord}{2};

            %%- create directories to the channel files per word pair
            wordoneDir = fullfile(sessionBlockDir, wordone);
            filesone = dir(fullfile(wordoneDir));
            filesone = {filesone(3:end).name};
            wordtwoDir = fullfile(sessionBlockDir, wordtwo);
            filestwo = dir(fullfile(wordtwoDir));
            filestwo = {filestwo(3:end).name};

            %%- loop through channels of data
            firstPairFeatureMat = [];
            secondPairFeatureMat = [];
            for iChan=4:4%length(filesone) % loop through and open up all channels
                % Load in the data struct for each word pair per channel
                fileOnePath = fullfile(wordoneDir, filesone{iChan});
                dataOne = load(fileOnePath);
                dataOne = dataOne.data;
                fileTwoPath = fullfile(wordtwoDir, filestwo{iChan});
                dataTwo = load(fileTwoPath);
                dataTwo = dataTwo.data;

                if DEBUG,
%                 fileOnePath
%                 wordone
%                 dataOne
%                 dataTwo
                end

                % concatenate all the freq. vectors that are already 500 ms
                % windowed and 100 ms overlap
                if isempty(firstPairFeatureMat),
                    firstPairFeatureMat = dataOne.powerMatZ;
                    secondPairFeatureMat = dataTwo.powerMatZ;
                else
                    firstPairFeatureMat = cat(2, firstPairFeatureMat, dataOne.powerMatZ);
                    secondPairFeatureMat = cat(2, secondPairFeatureMat, dataTwo.powerMatZ);
                end
                
                % plot avge spectrograms per channel for this specific wordpair
%                 figure
%                 imagesc(squeeze(mean(dataOne.powerMatZ,1)));
%                 hold on
%                 colormap('jet');
%                 title(['Same Pairs for block ', blocks{iBlock}])
%                 set(gca,'tickdir','out','YDir','normal');
%                 colorbar();
%                 
%                 figure
%                 imagesc(squeeze(mean(dataTwo.powerMatZ,1)));
%                 hold on
%                 colormap('jet');
%                 title(['Diff Pairs for block ', blocks{iBlock}])
%                 set(gca,'tickdir','out','YDir','normal');
%                 colorbar();
            end
            
            % create block spectrogram of all words
            if isempty(blockSpect)
                blockSpect = firstPairFeatureMat;
            else
                blockSpect = cat(1, blockSpect, firstPairFeatureMat);
            end
            
            % check the concatenated features
%             size(firstPairFeatureMat)
%             size(secondPairFeatureMat)
            
%             figure
%             imagesc(squeeze(mean(firstPairFeatureMat,1)));
%             hold on
%             colormap('jet');
%             title(['Same Pairs for block ', blocks{iBlock}])
%             set(gca,'tickdir','out','YDir','normal');
%             colorbar();
% 
%             figure
%             imagesc(squeeze(mean(secondPairFeatureMat,1)));
%             hold on
%             colormap('jet');
%             title(['Diff Pairs for block ', blocks{iBlock}])
%             set(gca,'tickdir','out','YDir','normal');
%             colorbar();
            
            %%- build events X features X time matrix with events being
            %%lined up to be compared
            for i=1:size(firstPairFeatureMat, 1),
                diffPairFeatureMat1 = cat(1, diffPairFeatureMat1, repmat(firstPairFeatureMat(i,:,:), size(secondPairFeatureMat, 1), 1, 1));
                diffPairFeatureMat2 = cat(1, diffPairFeatureMat2, secondPairFeatureMat(:, :, :));
            end
        end
        
        figure;
        imagesc(squeeze(mean(blockSpect, 1)));
        colormap('jet');
        colorbar();
        set(gca,'tickdir','out','YDir','normal');
        filesone{iChan}
        
        size(diffPairFeatureMat1)
        size(diffPairFeatureMat2)
    end
end