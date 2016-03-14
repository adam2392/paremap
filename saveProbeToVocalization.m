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
    waveletFreqs, trigType, thisChan, thisChanStr)
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
    
    buffer_powerMatZ = zeros(length(events), length(rangeFreqs));
    %%- Average From probe on -> vocalization [0: responseTime];
    responseTime = floor([events.responseTime]); % extract response Time's
    
    for i=1:length(events)
        buffer_powerMatZ(i,:) = mean(newPowerMatZ(i,:,2500:2500+responseTime(i)+1),3);
    end
    newPowerMatZ = buffer_powerMatZ;
    
    %%- plotting for each trigger
    uniqueTrigType = unique(trigType);          % get all unique triggers (e.g. all probe words)
            
    clear buffer_powerMatZ 
    %%- Save this new power matrix Z
    data.uniqueTrigType = uniqueTrigType; % store all the unique trigger types
    data.trigType = trigType;             % store the trigger type per event
    data.powerMatZ = newPowerMatZ;        % save the condensed power Mat
    data.chanNum = thisChan;           % store the corresponding channel number
    data.chanStr = thisChanStr;               % the string name of the channel
    data.freqBandYtick = freqBandYticks;
    data.freqBandYlabel = freqBandYtickLabels;
    
    filename = strcat(dataDir, num2str(thisChan), '_', thisChanStr, '_probeToVocal'); 
    save(filename, 'data');
end