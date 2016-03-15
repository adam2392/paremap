%%%%% Function: saveProbeToVocalization.m
%%%%% Description: Function to help save powerMatZ and meta data with a
%%%%% frequency binned (7 bins currently) and time binned from
%%%%% probeOn to Vocalization time.
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
function saveProbeToVocalization(events, powerMatZ, freqBandAr, ...
    waveletFreqs, waveT, trigType, thisChan, thisChanStr)
    %% Save Data As Frequency/ProbeToVocalization Binned
    dataDir = 'condensed_data/';
    
    % Create frequency band y ticks and ylabels
    freqBandYticks  = unique([freqBandAr(1:7).rangeF]);
    for iFB=1:length(freqBandYticks), 
        freqBandYtickLabels{iFB} = sprintf('%.0f Hz', freqBandYticks(iFB)); 
    end
    
    buffer_powerMatZ = squeeze(powerMatZ);
    
    %%- Call function to bin on freq. based on wavelet freqs. we have
    rangeFreqs = [freqBandAr.rangeF];
    rangeFreqs = reshape(rangeFreqs, 2, 7)';
    newPowerMatZ = freqBinSpectrogram(buffer_powerMatZ, rangeFreqs, waveletFreqs);
    
    % if no time domain needed (averaged across probeon -> vocalization
%     buffer_powerMatZ = zeros(length(events), length(rangeFreqs));
  
    
    %%- Average From probe on -> vocalization [0: responseTime];
    responseTime = floor([events.responseTime]); % extract response Time's
    
    WinLength = 100; % 100 ms
    Overlap = 50;    % overlap we want to increment
    NumWins = ceil(max(responseTime) / (WinLength-Overlap))+1;
    
%     buffer_powerMatZ(:,:) = mean(newPowerMatZ(:,:,2500:2500+responseTime(i)+1),3);
%     newPowerMatZ = buffer_powerMatZ;
    
    % condensed power mat on freq and time scale
    buffer_powerMatZ = newPowerMatZ(:,:,2500:2500+NumWins*(WinLength-Overlap)-1);
    
    % record the new response times, so we know when vocalization occurs in
    % condensed power matrix
    newResponseTimes = responseTime/Overlap;

    NumWins = size(buffer_powerMatZ,3) / (WinLength-Overlap) - 1;
    %%- Call function to bin on time based on winLength and Overlap and NumWins
    newPowerMatZ = timeBinSpectrogram(buffer_powerMatZ, NumWins, WinLength, Overlap);
     
    timeZero = find(waveT==0,1)/Overlap + 1; % index of timezero in bins
    
    % create time vector that is binned and still centered at 0
    binnedWaveT = 1:size(newPowerMatZ,3) - timeZero;
    
    clear buffer_powerMatZ 
    %%- Save this new power matrix Z
    data.trigType = trigType;             % store the trigger type per event
    data.powerMatZ = newPowerMatZ;        % save the condensed power Mat
    data.waveT = binnedWaveT;             % save the binned Wave T
    data.chanNum = thisChan;           % store the corresponding channel number
    data.chanStr = thisChanStr;               % the string name of the channel
    data.responseTime = newResponseTimes;%store the response time relative to condensed power matz
    data.freqBandYtick = freqBandYticks;
    data.freqBandYlabel = freqBandYtickLabels;
    
    filename = strcat(dataDir, num2str(thisChan), '_', thisChanStr, '_probeToVocal'); 
    save(filename, 'data');
end