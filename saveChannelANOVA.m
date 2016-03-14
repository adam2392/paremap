%%%%% Function: extractChannelFeature.m
%%%%% Description: Compute a distance metric from the condensed_data
%%%%% directory saved output. The data is in time binned and freq. binned
%%%%% matrix.
%%%%% 
%%%%% Input: 
%%%%% channel_num - the channel number we want to extract (e.g.
%%%%% 1-96)
%%%%% triggers - a list of triggers for all the events 
%%%%% powerMatZ - a power matrix, #eventsX#freqsXtime
%%%%% 
%%%%% Output:
%%%%% anovaMat - vector of datapoints (feature) for this channel
%%%%%
function saveChannelANOVA(powerMatZ, freqBandAr, trigType, ...
    thisChan, thisChanStr, WinLength, Overlap)
    %% RUN ANOVA 
    % squeeze spectrogram into eventsXfreqsXtime
    spectMat = squeeze(powerMatZ);
    
    %%%%%%%%%%%%%% Need to bin time and frequency first %%%%%%%%%%%%%%
    NumWins = size(spectMat,3) / (WinLength-Overlap) - 1;
    
    spectMat = timeBinSpectrogram(spectMat, NumWins, WinLength, Overlap);
    
    %%- Call function to bin on freq. based on 
    rangeFreqs = [freqBandAr.rangeF];
    rangeFreqs = reshape(rangeFreqs, 2, 7)';
    spectMat = freqBinSpectrogram(spectMat, rangeFreqs, waveletFreqs);
    
    size(spectMat)
    
    %% GET TRIGGERS WE WANT
    sampEventsMeta = events;  % includes assocaited + and *
    probeWords = {sampEventsMeta.probeWord};
    targetWords = {sampEventsMeta.targetWord};

    anovaPowMat = {}; % create a cell array to store all powerMatrices for certain trigger
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
                disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'BRICK PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
            case 'CLOCK'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'CLOCK PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'JUICE'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'JUICE PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'PANTS'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'PANTS PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
           case 'GLASS'
                tempevents = sampEventsMeta(strcmp(probeWords, THIS_TRIGGER));
                disp(['Looking at trigger: ', THIS_TRIGGER]);
                metaYstr = 'GLASS PROBE';

                tempInd = find(strcmp({sampEventsMeta.probeWord}, THIS_TRIGGER));
                %%- get all the unique targetwords for BRICK probeword
                targets = unique({tempevents.targetWord});
            otherwise
                error('no event trigger selected');
        end %end of switch
        %%- Create groups for ANOVA
        anovaPowMat{i} = spectMat(tempInd,:,:);
    end
    
    %% Actually Run ANOVA For This Channel
    anovaMat = zeros(size(spectMat,2), size(spectMat,3));
    for freq=1:size(spectMat,2)
        for time=1:size(spectMat,3)
            y = [];
            groups = [];
            for i=1:length(anovaPowMat)
                % create vector of events we want to test
                y = [y; anovaPowMat{i}(:,freq,time)];
                
                % set groups
                group = ones(size(anovaPowMat{i}(:,freq,time)))*i;
                groups = [groups; group];
            end
            
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
    
    % Create frequency band y ticks and ylabels
    freqBandYticks  = unique([freqBandAr(1:7).rangeF]);
    for iFB=1:length(freqBandYticks), 
        freqBandYtickLabels{iFB} = sprintf('%.0f Hz', freqBandYticks(iFB)); 
    end
    
    %%- Save this new power matrix Z
    anovaData.trigType = trigType;             % store the trigger type per event
    anovaData.anovaMat = anovaMat;        % save the condensed power Mat
    anovaData.chanNum = thisChan;           % store the corresponding channel number
    anovaData.chanStr = thisChanStr;               % the string name of the channel
    anovaData.freqBandYtick = freqBandYticks;
    anovaData.freqBandYlabel = freqBandYtickLabels;
    
    filename = strcat(dataDir, num2str(thisChan), '_', thisChanStr, '_anovaProbes'); 
    save(filename, 'anovaData'); 
end