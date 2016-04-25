WinLength = 100; % 100 ms
Overlap = 50;    % overlap we want to increment
NumWins = size(squeeze(powerMatZ),3) / (WinLength-Overlap) - 1;

%%- Call function to bin on time based on winLength and Overlap and NumWins
newPowerMatZ = timeBinSpectrogram(squeeze(powerMatZ), NumWins, WinLength, Overlap);
%%- Call function to bin on freq. based on wavelet freqs. we have
rangeFreqs = [freqBandAr.rangeF];
rangeFreqs = reshape(rangeFreqs, 2, 7)';
newPowerMatZ = freqBinSpectrogram(newPowerMatZ, rangeFreqs, waveletFreqs);

% data directory to save the data
dataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/blocks/';

sampEventsMeta = events;  % a buffer copy of the events struct
probeWords = {sampEventsMeta.probeWord};
targetWords = {sampEventsMeta.targetWord};

%%- LOOP through session, block, probe word and target word
%%- Get indices of all session blocks
sessions = unique({sampEventsMeta.sessionName});    % session0-2?
blocks = unique({sampEventsMeta.blocknumber});      % block0-5
for seshI=1:length(unique(sessions)), % loop thru session
    for blockI=1:length(unique(blocks)), % loop thru each block
        % get the event indices for this specific session block
%                 session_block_indices = find(strcmp({sampEventsMeta.sessionName}, sessions(seshI)) & ...
%                     strcmp({sampEventsMeta.blocknumber}, blocks(blockI)));
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
%                     tempevents = sampEventsMeta(strcmp({sampEventsMeta(session_block_indices).probeWord}, THIS_TRIGGER));
            targets = unique({tempevents.targetWord});

            %%- 02: GO THROUGH EACH TARGETWORD FOR THIS PROBEWORD (THISTRIGGER)
            %%-> Store all unique probe/target word pairs
            for j=1:length(targets) % loop through each unique trigger for a specific probeword
                % find event indices for this trigger matched with a specific
                % targetword
                targetWord = targets{j};

                % match probe, target and session and block
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

                clear data thisPowMat
            end
        end
    end
end