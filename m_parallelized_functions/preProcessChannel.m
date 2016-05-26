function preProcessChannel(iChan)    
    load('tempworkspace');
    powerMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    powerMatZ = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    phaseMat  = zeros(numChanPrealloc, length(eventTrigger), length(waveletFreqs), DurationMS);
    
    thisChan = chanList(iChan,:);   % the channel to use in this loop (e.g. 48)
    thisChanStr = chanStr{iChan};
    strStart    = sprintf('\n STEP 5.%d -- Grab %d/%d: %s', iChan, iChan, numChannels, thisChanStr );  strStart(end+1:35)=' '; %buffer length so everything lines up
    fprintf('%s', strStart);   tic;
    
    %%- gete_ms: get the eegWaveV
    % eegwaveform for each event over the duration of time for a certain channel
    eegWaveV = gete_ms(thisChan,eventTrigger,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,preFiltFreq,preFiltType,preFiltOrder,resampledrate);

    % notch filter to eliminate 60 Hz noise
    fprintf(' [%.1f sec] --> notch filt\n', toc); tic;
    eegWaveV = buttfilt(eegWaveV,[59.5 60.5],resampledrate,'stop',1); %-filter is overkill: order 1 --> 25 dB drop (removing 5-15dB peak)
    
    if ~ROBUST_SPEC % OPTION 1: perform wavelet spectral analysis
        %%- multiphasevec3: get the phase and power
        % power, phase matrices for events x frequency x duration of time for each channel
        fprintf(' [%.1f sec] --> freq decomp', toc); tic;
        [rawPhase,rawPow] = multiphasevec3(waveletFreqs,eegWaveV,resampledrate,waveletWidth);
        fprintf(' [%.1f sec] --> save', toc);  tic;
        fprintf('\n');

        %%- REMOVE LEADING/TRAILING buffer areas from power, phase,
        %%eegWave, timeVector
        rawPow   = rawPow(:,:,BufferMS+1:end-BufferMS);
        rawPhase = rawPhase(:,:,BufferMS+1:end-BufferMS);
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
        eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
        if length(eegWaveT)<size(eegWaveV,2), % error check on time vs. voltage length
            fprintf('wave time length off'); 
            eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
        end
        % x-axis of time series
        waveT = eegWaveT;

        % temp indicies
        iEv = 1:length(eventTrigger); % # of events
        iT  = 1:size(eegWaveV,2); % # of time points
        iF  = 1:length(waveletFreqs); % # of freqs.
        iChanSave = 1;

        % chan X event X freq X time
        % make power 10*log(power)
        powerMat(iChanSave,iEv,iF,iT) = 10*log10(rawPow);
        phaseMat(iChanSave,iEv,iF,iT) = rawPhase;

    %     for each eegfile stem, z-score each channel and frequency
        fprintf(' [%.1f sec] --> z-score', toc);  tic;
        stemList = unique({eventTrigger.eegfile});
        
        % indices of the powerMat to Z-score wrt
        for iStem=1:length(stemList),
            fprintf('.');
            iEvStem = find(strcmp({eventTrigger.eegfile}, stemList{iStem}));
            
            length(iEvStem)
            for iF = 1:length(waveletFreqs),
                allVal = reshape(squeeze(powerMat(iChanSave,iEvStem,iF,iT)),length(iEvStem)*length(iT),1); %allVal for particular chan and freq
                mu = mean(allVal); stdev = std(allVal);

                % create the power matrix
                powerMatZ(iChanSave,iEvStem,iF,iT) = (powerMat(iChanSave,iEvStem,iF,iT)-mu)/stdev;

                if sum(isnan(powerMatZ(iChanSave,iEvStem,iF,iT)))>0
                    keyboard;
                end
            end
        end
        % set two paramters from robust spectrotemp to 0 
        tWin = 0;
        freq = 0;
        
        fprintf(' [%.1f sec]', toc); tic;
%         clear rawPow rawPhase
        disp('powerMatZ, powerMat and phaseMat are created')
    else % OPTION 2: Perform robust spectrotemporal pursuit instead 
        %%- REMOVE LEADING/TRAILING BUFFER REGION FOR EEGWAVE
        eegWaveV = eegWaveV(:,BufferMS+1:end-BufferMS); % remove buffer area
        eegWaveT = (OffsetMS:DurationMS+OffsetMS-1)/1000; 
        if length(eegWaveT)<size(eegWaveV,2), % error check on time vs. voltage length
            fprintf('wave time length off'); 
            eegWaveT = (OffsetMS:DurationMS+OffsetMS)/1000;  
        end
        
        %%- ROBUST SPECTROTERMPORAL PURSUIT PARAMETERS
        % initialize arrays for storing 
        powerMat = zeros(length(eventTrigger), WINDOW/2, length(eegWaveV)/WINDOW); % #EVENTS X #FREQ X #TIME
        powerMatZ = zeros(length(eventTrigger), WINDOW/2, length(eegWaveV)/WINDOW);
        
        % generate robust spectrogram for each event
        tic;
        for iEv=1:length(eventTrigger),
            [xEst,freq,tWin,iter] = specPursuit(eegWaveV(iEv,:), FS, WINDOW, ALPHA);
            size(xEst)
            xEst = 20*log10(abs(xEst)); % convert to power spectrum
            xEst = reshape(xEst, 1, size(xEst,1), size(xEst,2));
            powerMat(iEv,:,:) = xEst;
        end
        fprintf(' [%.1f sec] --> robust spect pursuit', toc);
        
        % reset tWin var to reflect over seconds
        tWin = tWin - timeZero*1000;
        
        %%- Zscore all power for robust spectral pursuit
        fprintf(' [%.1f sec] --> z-score robust spec pursuit', toc); tic;
        for iEvent =1:size(powerMat,1),
            for iF = 1:size(powerMat,2),
                fixedVal = powerMat(iEvent, iF, 1:4); %allVal for particular chan and freq
                mu = mean(fixedVal); stdev = std(fixedVal);

                % create the power matrix
                powerMatZ(iEvent, iF, :) = (powerMat(iEvent, iF, :)-mu)/stdev;

                if sum(isnan(powerMatZ(iEvent,iF,:)))>0
                    keyboard;
                end
            end
        end
        fprintf(' [%.1f sec]', toc);
    end

    % create vector of the actual seconds in time axis for the powerMat
    % (since its time binned)...
    OVERLAP = 100;
    WINSIZE = 500;
    if tWin == 0, % if not set yet
        tWin = (LOWERTIME) :OVERLAP/FS: (UPPERTIME)-WINSIZE/FS;
    end
    timeZero = abs(0-(LOWERTIME))/(OVERLAP/FS);

  % remake powerMatZ to the points that we want (before probe on -> 3.5
    % seconds later
    powerMatZ = squeeze(powerMatZ); % only get the powerMatZ time points we want... (1001 - 2000+3500) = -1.0 seconds -> 3.5 seconds
%     powerMatZ = powerMatZ(:,:,LOWERTIME:UPPERTIME);
    
    if DEBUG,
        size(powerMatZ)
    end
    %% TIME BIN POWERMATZ WITH WINDOWSIZE AND OVERLAP
    addpath('../m_oldAnalysis_anovaANDsinglechannel/');
    WINDOWSIZE = 500; % in milliseconds
    OVERLAP = 100;    % in milliseconds
    powerMatZ = timeBinSpectrogram(powerMatZ, WINDOWSIZE, OVERLAP);
 
    if DEBUG,
        size(powerMatZ)
    end
    
    %% FREQUENCY BIN WITH FREQUENCY BANDS
    rangeFreqs = reshape([freqBandAr.rangeF], 2, 7)';
    waveletFreqs = waveletFreqs;
    powerMatZ = freqBinSpectrogram(powerMatZ, rangeFreqs, waveletFreqs);

    if DEBUG,
        size(powerMatZ)
    end
    % SPLIT INTO SESSIONS AND BLOCKS
    subjSessions = unique({events.sessionName}); % e.g. sessions 0, 1, 2
    subjBlocks = unique({events.blocknumber});   % e.g. blocks 0,1,2,3,4,5
    
    for iSesh=1:length(subjSessions),
        for iBlock=1:length(subjBlocks),
            sessionBlockIndices = strcmp({events.sessionName}, subjSessions(iSesh)) & ...
                                    strcmp({events.blocknumber}, subjBlocks(iBlock));
                                
            %%- LOOP THROUGH PROBEWORDS
            probeWords = unique({events(sessionBlockIndices).probeWord});
            for iProbe=1:length(probeWords),
                THIS_PROBE = probeWords{iProbe};
                
                % events for this probeWord and their targetWords
                probeIndices = strcmp({events.probeWord}, THIS_PROBE);
                tempEvents = events(probeIndices & sessionBlockIndices);
                targetWords = unique({tempEvents.targetWord});
                
                %%- LOOP THROUGH TARGETWORDS FOR EACH PROBEWORD
                for iTarget=1:length(targetWords),
                    THIS_TARGET = targetWords{iTarget};
                    
                    % match probe, target, session and block
                    eventIndices = find(strcmp({events.probeWord}, THIS_PROBE) & ...
                                    strcmp({events.targetWord}, THIS_TARGET) & ...
                                    sessionBlockIndices);
                    sessionBlockWordPairEvents = events(eventIndices);
                    
                    % vars to save
                    blockNum = unique({sessionBlockWordPairEvents.blocknumber});
                    sessionNumber = sessionBlockWordPairEvents(1).sessionNum;
                    thisPowMat = powerMatZ(eventIndices,:,:);
                    eegWaveToSave = eegWaveV(eventIndices, :);
                    freqBandYtick = 1:length(freqBandYticks);
                    freqBandYlabel = {freqBandAr.name};
                    vocalization = ceil([sessionBlockWordPairEvents.responseTime]/OVERLAP);
                    %% SAVE PROCESSED DATA IN A MATLAB STRUCT
                    if SAVE,
                        data.probeWords = THIS_PROBE;                   % the probe words for all events in this struct
                        data.targetWords = THIS_TARGET;                 % the target words for all events in this struct
                        data.sessionNum = sessionNumber;                % the session number
                        data.blockNum = blockNum;                       % the block number
                        data.eegWaveV = eegWaveToSave;       % eeg wave form
                        data.eegWaveT = eegWaveT;                       % time series for eeg voltage wave
                        data.chanNum = thisChan;                        % store the corresponding channel number
                        data.chanStr = thisChanStr;                     % the string name of the channel
                        data.freqBandYtick = freqBandYtick;            % store frequency bands if using wavelet transform
                        data.freqBandYlabel = freqBandYlabel;
                        data.descriptor = 'Initial processing -1 seconds to 5 seconds after probeWordOn. Time binned with 500ms window and 100ms overlap';
                        data.timeZero = timeZero; %ceil((TIMEZERO-LOWERTIME)/OVERLAP);
                        data.vocalization = data.timeZero + vocalization;
                        data.powerMatZ = thisPowMat;            % save the condensed power Mat Z-scored
                        data.waveT = tWin;                      % ROBUSTSPECT: save the binned Wave T
                        data.freq = freq;                       % ROBUSTSPECT: save the frequency points

                        %%- SAVING DIR PARAMETERS
                        if ROBUST_SPEC,
                            TYPE_SPECT = 'robust_spec';
                        else
                            TYPE_SPECT = 'morlet_spec';
                        end
                        if VOCALIZATION
                            TYPE_SPECT = strcat(TYPE_SPECT, '_vocalization');
                        elseif MATCHWORD
                            TYPE_SPECT = strcat(TYPE_SPECT, '_matchword');
                        end

                        chanFileName = strcat(num2str(thisChan), '_', thisChanStr, '_', TYPE_SPECT);
                        wordpair_name = strcat(THIS_PROBE, '_', THIS_TARGET);

                        % data directories to save data into
                        workDir = '/Users/liaj/Documents/MATLAB/paremap';
                        homeDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/';
                        homeDir = '/home/adamli/paremap/';
                        dataDir = strcat('condensed_data_', subj);
                        typeTransformDir = fullfile(homeDir, dataDir, TYPE_SPECT);
                        fileDir = fullfile(typeTransformDir, subjSessions{iSesh}, subjBlocks{iBlock}, wordpair_name);
                        chanFilePath = fullfile(fileDir, chanFileName);; 

                        if ~exist(fileDir)
                            mkdir(fileDir);
                        end
                        save(chanFilePath, 'data');  
                    end
                end % loop through target
            end % loop through probe
        end % loop through block 
    end % loop through session
end