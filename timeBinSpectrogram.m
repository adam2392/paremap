%%% function timeBinSpectrogram
%%% Input: 
%%%  - spectMat:     #events X #freqs. X #timepoints = 3D array
%%%  - numTimeWins:  The number of time windows we want
%%% Output:
%%%  - newSpect: the new spectrogram with (#events X #freqs X numTimeWins)
%%%  size array
function newSpect = timeBinSpectrogram(spectMat, numTimeWins, WinLength, Overlap)

    % initialize new power matrix
    newSpect = zeros(size(spectMat,1),size(spectMat,2),numTimeWins); 

    % loop through the number of events 
    for i=1:size(spectMat,1),
        clear avgpowermat % clear so that each timewindow, creates new var
        
        for j=1:numTimeWins % loop through and bin the time domain into windows
            % index through the time domain of the Z power matrix
            index = 1+(j-1)*Overlap;

            % create buffer variable for power for each time window
            if j == numTimeWins
                eventpowerMat = spectMat(i,:,index:end);
            else
                eventpowerMat = spectMat(i,:,index:index+WinLength);
            end

            % average the power in time windows and append to vector
            if ~exist('avgpowermat')
                avgpowermat = mean(eventpowerMat,3);
            else
                avgpowermat = [avgpowermat; mean(eventpowerMat,3)];
            end
        end
        % append avgpowermat to new matrix
        newSpect(i,:,:) = avgpowermat';
    end
    
    % check if return array dimensions are correct
    if ~isequal(size(newSpect),[size(spectMat,1),size(spectMat,2),numTimeWins]),
        disp('error in timeBinSpectrogram.m')
    end
end