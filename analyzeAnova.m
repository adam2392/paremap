%%%%% Function: analyzeANOVA.m
%%%%% Description: Saves anova p-value matrix using precomputed data from
%%%%% 100ms time binned data. This saves time. Then this correspondingly
%%%%% saves a struct that contains the anova output that can be analyzed
%%%%% using analyzeAnovaMat.m.
%%%%% 
%%%%% Input: 
%%%%% data directory saved from previous running of script at 100ms/50ms
%%%%% overlap
%%%%% 
%%%%% Output:
%%%%% saves data as 'filename' variable set in script
%%%%%
anovaDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/freq_probeToVocal_100msbinned/';
ext = '*.mat';
files = dir(strcat(anovaDir, ext));
files = {files.name};
file = strcat(anovaDir, files{1});
data = load(file);
data = data.data;
%%-array of frequency bands
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

freqBandYticks = data.freqBandYtick;
freqBandYlabels = data.freqBandYlabel;
thisChan = data.chanNum;
thisChanStr = data.chanStr;

subj = 'NIH034';

%%- event data 
eegRootDirWork = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg/';     % work
eegRootDirHome = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';  % home

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

% loop through each channel file
for filei=1:length(files)
    
    file = strcat(anovaDir, files{filei});
    data = load(file);
    data = data.data;
    
    disp(['on file: ', file])
    %% EXTRACT DAT THAT WE WANT
    thisChan = data.chanNum;
    thisChanStr = data.chanStr;
    powerMatZ = data.powerMatZ;
    trigType = data.trigType;
    responseTimes = data.responseTime;
    
    %% RUN ANOVA 
    %%- Data is already time and frequency binned
    % squeeze spectrogram into eventsXfreqsXtime if necessary
    spectMat = squeeze(powerMatZ);
    size(spectMat);
    
    %% GET TRIGGERS WE WANT
    sampEventsMeta = events;  % includes assocaited + and *
    probeWords = {sampEventsMeta.probeWord};
    targetWords = {sampEventsMeta.targetWord};
    TRIGGER_TYPES   = {'BRICK','CLOCK','GLASS','JUICE','PANTS'}; %-fixation and blockStart not ready for physio analysis

    anovaPowMat = {}; % create a cell array to store all powerMatrices for certain trigger
    anovaGroups = {}; % create cell array to hold the groups string
    index = 0;
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
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'BRICK PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
            case 'CLOCK'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'CLOCK PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'JUICE'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'JUICE PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'PANTS'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'PANTS PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'GLASS'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
%                 disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'GLASS PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
            otherwise
                error('no event trigger selected');
        end %end of switch
        
        %%- 02: GO THROUGH EACH TARGETWORD FOR THIS PROBEWORD (THISTRIGGER)
        for j=1:length(targets) % loop through each unique trigger for a specific probeword
            % find event indices for this trigger matched with a specific
            % targetword
            targetWord = targets{j};
            eventInd = find(strcmp({sampEventsMeta.probeWord},THIS_TRIGGER) & strcmp({sampEventsMeta.targetWord},targetWord));
            
            index = index + 1;
            %%- Create groups for ANOVA
            anovaGroups{index} = strcat(THIS_TRIGGER, '_', targetWord);
            anovaPowMat{index} = spectMat(eventInd,:,:);
        end
        
        %%- Create groups for ANOVA
%         anovaPowMat{i} = spectMat(tempInd,:,:);
    end
    
    %%- RUN THESE 2 LINES TO DETERMINE WHICH INDICES YOU WANT IN 'Y'
%     anovaPowMat
%     anovaGroups
    first = 3;
    second = 7;
    %% Actually Run ANOVA For This Channel
    anovaMat = zeros(size(spectMat,2), size(spectMat,3));
    for freq=1:size(spectMat,2)
        for time=1:size(spectMat,3)
            y = [];
            groups = [];
            
            y = [anovaPowMat{first}(:,freq,time); anovaPowMat{second}(:,freq,time)];
            group = ones(size(anovaPowMat{first}(:,freq,time)));
            groups = [ones(size(anovaPowMat{first}(:,freq,time)));...
                      ones(size(anovaPowMat{second}(:,freq,time)))*2];
            
            %%- loop through every group
%             for i=1:length(anovaPowMat)
%                 % create vector of events we want to test
%                 y = [y; anovaPowMat{i}(:,freq,time)];
%                 
%                 % set groups
%                 group = ones(size(anovaPowMat{i}(:,freq,time)))*i;
%                 groups = [groups; group];
%             end
            
            % compute p-value for ANOVA
            p = anovan(y, groups, 'display','off');
            %%%% 2nd group with target words,...
            
            % add to ANOVA matrix of p-values
            anovaMat(freq,time) = p;
        end
    end
%     anovaMat(anovaMat > 0.05) = 1;

    %% Save Data As Frequency/ProbeToVocalization Binned
    dataDir = 'condensed_data/';
    clear data
    
    %%- Save this new power matrix Z
    data.trigType = trigType;             % store the trigger type per event
    data.anovaMat = anovaMat;        % save the condensed p-value matrix
    data.chanNum = thisChan;           % store the corresponding channel number
    data.chanStr = thisChanStr;               % the string name of the channel
    data.freqBandYticks = freqBandYticks;
    data.freqBandYlabels = freqBandYlabels;
    data.responseTimes = responseTimes;
    
    filename = strcat(dataDir, num2str(thisChan), '_', thisChanStr, ...
        '_anovaProbeOnToVocalization_', anovaGroups{first}, anovaGroups{second}); 
    save(filename, 'data'); 
    
    clear data
end