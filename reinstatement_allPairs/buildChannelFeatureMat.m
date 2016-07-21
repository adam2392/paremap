function [wordPairFeatureMat1, wordPairFeatureMat2] = buildChannelFeatureMat(channels, wordPairDir1, arg3)
    % inputs: channels, wordPairDir
    wordPairFeatureMat1 = [];
    wordPairFeatureMat2 = [];
    
    if nargin == 3
        wordPairDir2 = arg3;
    end
    
    %%- either build up feature matrix for 1 wordPair or 2
    if ~exist('wordPairDir2', 'var')
        for iChan=1:length(channels)
            data = load(fullfile(wordPairDir1, channels{iChan}));
            data = data.data;

            if isempty(wordPairFeatureMat1),
                wordPairFeatureMat1 = data.powerMatZ;
            else
                wordPairFeatureMat1 = cat(2, wordPairFeatureMat1, data.powerMatZ);
            end
        end
    else
       for iChan=1:length(channels) % loop through and open up all channels
            % Load in the data struct for each word pair per channel
            dataOne = load(fullfile(wordPairDir1, channels{iChan}));
            dataOne = dataOne.data;
            dataTwo = load(fullfile(wordPairDir2, channels{iChan}));
            dataTwo = dataTwo.data;

            % concatenate all the freq. vectors that are already 500 ms
            % windowed and 100 ms overlap
            if isempty(wordPairFeatureMat1),
                wordPairFeatureMat1 = dataOne.powerMatZ;
                wordPairFeatureMat2 = dataTwo.powerMatZ;
            else
                wordPairFeatureMat1 = cat(2, wordPairFeatureMat1, dataOne.powerMatZ);
                wordPairFeatureMat2 = cat(2, wordPairFeatureMat2, dataTwo.powerMatZ);
            end
        end 
    end
end