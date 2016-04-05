%%%% USED TO SAVE ROBUST SPECTROGRAM INTO SEPARATE BLOCKS/SESSIONS

clear
clc
close all

%% SUBJECT AND BLOCK SELECTION
subj = 'NIH034'; 
sessNum = [0:2];

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


eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end

% Either go through all the sessions, or a specific session
if sessNum == -1 | length(sessNum)>1, % all sessions
    disp('STEP 1: Going through all sessions')
    session = 'Meta Session [all]';
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/');
    sessStr = '[all]';
else                                  % one session
    disp('STEP 1: Going through one session')
    session = sprintf('session_%d',sessNum);
    behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/', session);
    sessStr = sprintf('[%d]',sessNum);
end

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

% only get the correct events
correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);
clear correctIndices 

%%- Load Robust Spect Results
robustDir = strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/robust_spec/');

ext = '*.mat';
files = dir(strcat(robustDir, ext));
files = {files.name};

for fileI=1:length(files),
    filename = strcat(robustDir, files{fileI});
    data = load(filename);
    data = data.data;
    
    % already freq/time binned
    newPowerMatZ = 10*log10(abs(data.powerMatZ));
    
    if sum(sum(sum(isnan(newPowerMatZ))))>0 || sum(sum(sum(isinf(newPowerMatZ))))>0
        filename
    end
    
    thisChan = data.chanNum;
    thisChanStr = data.chanStr;
    Overlap = 50;
        
    %%- Save
    %%- Call function to bin on freq. based on wavelet freqs. we have
    rangeFreqs = [freqBandAr.rangeF];
    rangeFreqs = reshape(rangeFreqs, 2, 7)';

    % data directory to save the data
    dataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/robustspec_blocks/';

    sampEventsMeta = events;  % a buffer copy of the events struct
    probeWords = {sampEventsMeta.probeWord};
    targetWords = {sampEventsMeta.targetWord};

    %%- LOOP through session, block, probe word and target word
    %%- Get indices of all session blocks
    sessions = unique({sampEventsMeta.sessionName});    % session0-2?
    blocks = unique({sampEventsMeta.blocknumber});      % block0-5
    for seshI=1:length(unique(sessions)),
        for blockI=1:length(unique(blocks)),
            disp(['Saving session ', seshI, ' and block ', blockI])
            % get the event indices for this specific session block
            session_block_indices = find(strcmp({sampEventsMeta.sessionName}, sessions(seshI)) & ...
                strcmp({sampEventsMeta.blocknumber}, blocks(blockI)));
            session_block_indices = strcmp({sampEventsMeta.sessionName}, sessions(seshI)) & ...
                strcmp({sampEventsMeta.blocknumber}, blocks(blockI));

            %%- Loop through only probewords for this sessionblock
            probeWords = unique({sampEventsMeta(session_block_indices).probeWord});
            for i=1:length(probeWords)
                THIS_TRIGGER = probeWords{i};

                % get the events for this probe word and get the
                % corresponding targets
                probeInd = strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER);
                tempevents = sampEventsMeta(probeInd & session_block_indices);
                targets = unique({tempevents.targetWord});
                
                %%- 02: GO THROUGH EACH TARGETWORD FOR THIS PROBEWORD (THISTRIGGER)
                %%-> Store all unique probe/target word pairs
                for j=1:length(targets) % loop through each unique trigger for a specific probeword
                    % find event indices for this trigger matched with a specific
                    % targetword
                    targetWord = targets{j};
                    eventInd = find(strcmp({sampEventsMeta.probeWord},THIS_TRIGGER) & ...
                        strcmp({sampEventsMeta.targetWord},targetWord) & session_block_indices);
                    metaEvents = events(eventInd);

                    %%- store each relevant power matrix
                    thisPowMat = newPowerMatZ(eventInd,:,:);

                    data.metaEvents = metaEvents;
                    data.powerMatZ = thisPowMat;
                    data.chanNum = thisChan;
                    data.chanStr = thisChanStr;
                    data.probeWord = THIS_TRIGGER;
                    data.targetWord = targetWord;
                    data.timeZero = 45; %%%%% ** MAGIC NUMBER BECAUSE 2.25-5.25
                    data.vocalization = data.timeZero + round([metaEvents.responseTime]/Overlap);

                    %%- save into this dir
                    wordpair_name = strcat(THIS_TRIGGER, '_', targetWord);
                    filename = strcat(dataDir, sessions{seshI}, '/', blocks{blockI}, '/', wordpair_name, '/', num2str(thisChan), '_', thisChanStr, '_groupSessionBlockData');

                    filedir = strcat(dataDir, sessions{seshI}, '/', blocks{blockI}, '/', wordpair_name, '/');
                    if ~exist(filedir)
                        mkdir(filedir);
                    end

                    save(filename, 'data');
                end
            end
        end
    end
end