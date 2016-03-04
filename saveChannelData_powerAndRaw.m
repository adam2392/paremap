function success = saveChannelData_powerAndRaw(data, rawdata, channelNum, chanString)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from paRemap %%%%
%
%
%   extraction designed for paRemap v4.2, which is implemented in pyEPL and was used for NIH030 and beyond
%           (earlier versions were implemented in psychoPy) and were used for NIH028 and 029... those session logs will need some tweaking to use with this extraction
%
%   training section NOT saved to session log... for earlier versions this MUST BE DELETED from the session log
%
%
%
%%%%%   create an event for every presented word and every response (no words from the training sections should be included)

% %%- TIME BIN Before saving
%     WinLength = 100; % 100 ms
%     Overlap = 50;
%     NumWins = size(powerMatZ,3) / (WinLength-Overlap) - 1;
% 
%     timeZero = find(waveT==0,1)/Overlap + 1; % index of timezero in bins
% 
%     % initialize new power matrix
%     newPowerMatZ = zeros(size(powerMatZ,1),size(powerMatZ,2),NumWins); 
% 
%     % window by time
%     for i=1:size(powerMatZ,1),
%         clear avgpowermat
%         % loop through and bin the time domain into windows
%         for j=1:NumWins %1:49:size(powerMatZ,3)-1
%             % index through the time domain of the Z power matrix
%             index = 1+(j-1)*50;
% 
%             if j == NumWins
%                 eventpowerMat = powerMatZ(i,:,index:end);
%             else
%                 eventpowerMat = powerMatZ(i,:,index:index+100);
%             end
% 
%             if ~exist('avgpowermat')
%                 avgpowermat = mean(eventpowerMat,3);
%             else
%                 avgpowermat = [avgpowermat; mean(eventpowerMat,3)];
%             end
%         end
%         % append avgpowermat to new matrix
%         newPowerMatZ(i,:,:) = avgpowermat';
%     end
% 
%     % create time vector that is binned and still centered at 0
%     binnedWaveT = 1:size(newPowerMatZ,3) - timeZero;
%     
%     data.powerMatZ = newPowerMatZ;
%     data.waveT = binnedWaveT;
    
    filename = strcat(num2str(channelNum), '_', chanString);
    save(filename, 'data', '-v7.3');
    
    filename = strcat(num2str(channelNum), '_', chanString, '_rawdata');
    save(filename, 'rawdata', '-v7.3');
    
    success = 1;
end