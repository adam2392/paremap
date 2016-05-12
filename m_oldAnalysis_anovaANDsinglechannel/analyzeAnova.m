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
dataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/freq_probeToVocal_100msbinned/';
% dataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/robust_spec/';

ext = '*.mat';
files = dir(strcat(dataDir, ext));
files = {files.name};
file = strcat(dataDir, files{1});
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
responseTime = data.responseTime;

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

% firstgroups = [1, 1, 2, 3];
% secondgroups = [7, 8, 8, 7];
% firstgroups = [7, 8, 3, 2];
% secondgroups = [10, 12, 11, 9];
firstgroups=[1, 2, 3, 4];  %same target words
secondgroups=[6, 7, 8, 9];

for groupind = 1:length(firstgroups)
   first = firstgroups(groupind)
   second = secondgroups(groupind)
    % loop through each channel file
    for filei=1:length(files)

        file = strcat(dataDir, files{filei});
        data = load(file);
        data = data.data;

        disp(['on file: ', file])
        %% EXTRACT DAT THAT WE WANT
        thisChan = data.chanNum;
        thisChanStr = data.chanStr;
        powerMatZ = data.powerMatZ;
        trigType = data.trigType;
    %     responseTimes = data.responseTime;

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
%         anovaPowMat
%         anovaGroups
    %     first = 1; %1/7, 1/8, 2/8, 3/7 (probed brick diff words completely)
    %     second = 7;
    % 7, 8, 3, 2 % flip flopped words
    % 10, 12, 11, 9
    % 
    % 1, 2, 3, 4  %same target words
    % 6, 7, 8, 9
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

                % compute p-value for ANOVA
                p = anovan(y, groups, 'display','off');
                %%%% 2nd group with target words,...

                % add to ANOVA matrix of p-values
                anovaMat(freq,time) = p;
            end
        end
    %     anovaMat(anovaMat > 0.05) = 1;

        %% Save Data As Frequency/ProbeToVocalization Binned
        newdataDir = strcat('condensed_data/anova/',...
            lower(strjoin(strsplit(anovaGroups{first}, '_'), '')), '_',  lower(strjoin(strsplit(anovaGroups{second}, '_'), '')), '/');
        if ~exist(newdataDir, 'dir')
            mkdir(newdataDir)
        end

        %%- Save this new power matrix Z
        newdata.trigType = trigType;             % store the trigger type per event
        newdata.anovaMat = anovaMat;        % save the condensed p-value matrix
        newdata.chanNum = thisChan;           % store the corresponding channel number
        newdata.chanStr = thisChanStr;               % the string name of the channel
        newdata.freqBandYticks = freqBandYticks;
        newdata.freqBandYlabels = freqBandYlabels;
        newdata.responseTimes = responseTime;
        newdata.firstgroup = anovaGroups{first};
        newdata.secondgroup = anovaGroups{second};

        filename = strcat(newdataDir, num2str(thisChan), '_', thisChanStr, ...
            '_anovaProbeOnToVocalization_', anovaGroups{first}, anovaGroups{second}); 
        save(filename, 'newdata'); 

        clear data
    end
end