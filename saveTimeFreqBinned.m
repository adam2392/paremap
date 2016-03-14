%%%%% Function: saveTimeFreqBinned.m
%%%%% Description: Function to help save powerMatZ and meta data with a
%%%%% frequency binned (7 bins currently) and time binned with winLength
%%%%% 
%%%%% Dependencies: timeBinSpectrogram and freqBinSpectrogram
%%%%% 
%%%%% Input: 
%%%%% channel_num - the channel number we want to extract (e.g.
%%%%% 1-96)
%%%%% triggers - a list of triggers for all the events 
%%%%% powerMatZ - a power matrix, #eventsX#freqsXtime
%%%%% freqBandAr - freq. Band Array
%%%%% waveletFreqs - an array of wavelet freqs.
%%%%% 
%%%%% Output:
%%%%% anovaMat - vector of datapoints (feature) for this channel
%%%%%
function saveTimeFreqBinned(powerMatZ, freqBandAr, waveletFreqs, ...
    trigType, thisChan, thisChanStr, WinLength, Overlap, waveT)

%%- Here is where I save the data
    %%- i) condense data in freq/time domain and ii) save it
    dataDir = 'condensed_data/';
    if ~exist(dataDir), mkdir(dataDir); end
    
    % Create frequency band y ticks and ylabels
    freqBandYticks  = unique([freqBandAr(1:7).rangeF]);
    for iFB=1:length(freqBandYticks), 
        freqBandYtickLabels{iFB} = sprintf('%.0f Hz', freqBandYticks(iFB)); 
    end
    
    % squeeze out the channel dimension of powerMatZ and time bin
    buffer_powerMatZ = squeeze(powerMatZ);
    %%- TIME BIN Before saving
    NumWins = size(buffer_powerMatZ,3) / (WinLength-Overlap) - 1;
    
    %%- Call function to bin on time based on winLength and Overlap and NumWins
    newPowerMatZ = timeBinSpectrogram(buffer_powerMatZ, NumWins, WinLength, Overlap);

    %%- Call function to bin on freq. based on 
    rangeFreqs = [freqBandAr.rangeF];
    rangeFreqs = reshape(rangeFreqs, 2, 7)';
    newPowerMatZ = freqBinSpectrogram(newPowerMatZ, rangeFreqs, waveletFreqs);
    
    timeZero = find(waveT==0,1)/Overlap + 1; % index of timezero in bins
    
    % create time vector that is binned and still centered at 0
    binnedWaveT = 1:size(newPowerMatZ,3) - timeZero;
    
    %%- plotting for each trigger
    uniqueTrigType = unique(trigType);          % get all unique triggers (e.g. all probe words)
    numUniqueTrig  = length(uniqueTrigType);    % get length of unique triggers
        
    %% Create power mat data struct
    %save data that can plot evoked and spectrogram
    data.uniqueTrigType = uniqueTrigType; % store all the unique trigger types
    data.trigType = trigType;             % store the trigger type per event
    data.powerMatZ = newPowerMatZ;        % save the condensed power Mat
    data.waveT = binnedWaveT;             % save the binned Wave T
    data.waveletFreqs = waveletFreqs;     % wavelet frequencies used
    data.chanNum = thisChan;           % store the corresponding channel number
    data.chanStr = thisChanStr;               % the string name of the channel
    data.freqBandYtick = freqBandYticks;
    data.freqBandYlabel = freqBandYtickLabels;
    data.timeZero = timeZero;
    %%%%% ADD INDICES OF EVENTS?
    
    filename = strcat(dataDir, num2str(thisChan), '_', thisChanStr); 
    save(filename, 'data');