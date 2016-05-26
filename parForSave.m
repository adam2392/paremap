function parForSave(chanFilePath, THIS_PROBE, THIS_TARGET, sessionNumber, blockNum, ...
    eegWaveToSave, eegWaveT, thisChan, thisChanStr, freqBandYtick, freqBandYlabel, timeZero, ...
    vocalization, thisPowMat, tWin, freq)
%%- Save this new power matrix Z-scored into data .mat file
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

    save(chanFilePath, 'data');  
end