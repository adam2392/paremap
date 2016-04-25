dataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/freq_probeToVocal_100msbinned/';
ext = '*.mat';
files = dir(strcat(dataDir, ext));
files = {files.name};
file = strcat(dataDir, files{1});
data = load(file);
data = data.data;

subj = 'NIH034';

%%- event data 
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/';  % home

% Determine which directory we're working with automatically
if     length(dir(eegRootDirWork))>0, eegRootDir = eegRootDirWork;
elseif length(dir(eegRootDirHome))>0, eegRootDir = eegRootDirHome;
else   error('Neither Work nor Home EEG directories exist! Exiting'); end
behDir=fullfileEEG(eegRootDir, subj, '/behavioral/paremap/');

%%-Load in the Events For This Task/Patient/Session
events = struct([]);                    %%- in functional form this is required so there is no confusion about events the function and events the variable
load(sprintf('%s/events.mat',behDir));  %%- load the events file
fprintf('Loaded %d events from %s\n', length(events), behDir);

% only get the correct events
correctIndices = find([events.isCorrect]==1);
events = events(correctIndices);
clear correctIndices 

newPowerMatZ = data.powerMatZ;

        % data directory to save the data
        dDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/groups/';
        
        sampEventsMeta = events;  % includes assocaited + and *
        probeWords = {sampEventsMeta.probeWord};
        targetWords = {sampEventsMeta.targetWord};
        %%- Loop through each probeword
        for i=1:length(TRIGGER_TYPES)
            THIS_TRIGGER = TRIGGER_TYPES{i}; % set the current probeword

            %%- 01: GET TRIGGER INDICES WE WANT
            switch THIS_TRIGGER,
                %%- For each probeword:
                % - find events with that probeword
                % - get the unique targetwords for that event
                case 'BRICK'
                    tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));            
                    tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                    %%- get all the unique targetwords for BRICK probeword
                    targets = unique({tempevents.targetWord});
                case 'CLOCK'
                    tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                    tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                    %%- get all the unique targetwords for BRICK probeword
                    targets = unique({tempevents.targetWord});
               case 'JUICE'
                    tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                    tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                    %%- get all the unique targetwords for BRICK probeword
                    targets = unique({tempevents.targetWord});
               case 'PANTS'
                    tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                    tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                    %%- get all the unique targetwords for BRICK probeword
                    targets = unique({tempevents.targetWord});
               case 'GLASS'
                    tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                    tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                    %%- get all the unique targetwords for BRICK probeword
                    targets = unique({tempevents.targetWord});
                otherwise
                    error('no event trigger selected');
            end

            %%- 02: GO THROUGH EACH TARGETWORD FOR THIS PROBEWORD (THISTRIGGER)
            %%-> Store all unique probe/target word pairs
            for j=1:length(targets) % loop through each unique trigger for a specific probeword
                % find event indices for this trigger matched with a specific
                % targetword
                targetWord = targets{j};
                eventInd = find(strcmp({sampEventsMeta.probeWord},THIS_TRIGGER) & strcmp({sampEventsMeta.targetWord},targetWord));
                metaEvents = events(eventInd);
                 
                %%- store each relevant power matrix
                thisPowMat = newPowerMatZ(eventInd,:,:);
                
                data.powerMatZ = thisPowMat;
                data.chanNum = thisChan;
                data.chanStr = thisChanStr;
                data.probeWord = THIS_TRIGGER;
                data.targetWord = targetWord;
                data.timeZero = 45; %%%%% ** MAGIC NUMBER BECAUSE 2.25-5.25
                data.vocalization = data.timeZero + round([metaEvents.responseTime]/Overlap);
                
                %%- save into this dir
                wordpair_name = strcat(THIS_TRIGGER, '_', targetWord);
                filename = strcat(dDir, wordpair_name, '/', num2str(thisChan), '_', thisChanStr, '_groupData');
                
                filedir = strcat(dDir, wordpair_name, '/');
                if ~exist(filedir)
                    mkdir(filedir);
                end
                
                save(filename, 'data');
            end
        end