%%% function timeBinSpectrogram
%%% Input: 
%%%  - spectMat:     #events X #freqs. X #timepoints = 3D array
%%%  - rangeFreqs:  A nx2 matrix that holds the ranges for each freq.
%%%  window
%%%  - waveletFreqs: The vector of wavelet frequencies for each of the
%%%  bands in spectMat
%%%
%%% Output:
%%%  - newSpect: the new spectrogram with (#events X numFreqWins X #timepoints)
%%%  size array
function newSpect = freqBinSpectrogram(spectMat, rangeFreqs, waveletFreqs)
    % get number of frequency windows we want
    numFreqWins = size(rangeFreqs, 1);

    % initialize new power matrix
    newSpect = zeros(size(spectMat,1),numFreqWins,size(spectMat,3)); 
    
    clear avgpowermat % clear so that each timewindow, creates new var

    %%- loop through numFreqWins rows
    for j=1:numFreqWins % loop through and bin the freq domain into windows
        lowerFreq = rangeFreqs(j, 1);
        if lowerFreq == 2
            lowerFreq = 0;
        end
        upperFreq = rangeFreqs(j, 2);

%         lowerFreq
%         upperFreq
        
        %%- go through indices in waveletFreqs and average those
        %%between lower and upper freq. -> append to eventpowerMat
        lowerInd = waveletFreqs >= lowerFreq;
        upperInd = waveletFreqs <= upperFreq;
        indices = lowerInd == upperInd; % binary selector index vector

        % create buffer variable for power for each time window
        eventpowerMat = spectMat(:, indices,:);

        % average the power in freq. windows and append to vector
        if ~exist('avgpowermat')
            avgpowermat = mean(eventpowerMat,2);
        else
            avgpowermat = cat(2, avgpowermat, mean(eventpowerMat,2));
        end
    end
%         size(avgpowermat)
    % append avgpowermat to new matrix
    newSpect(:,:,:) = avgpowermat;

    % check if return array dimensions are correct
    if ~isequal(size(newSpect),[size(spectMat,1),numFreqWins, size(spectMat,3),]),
        disp('error in timeBinSpectrogram.m')
    end
end