function [samePairFeatureMat1, samePairFeatureMat2] = createSamePairReinInput(channels, wordPairDir)
    samePairFeatureMat1 = [];
    samePairFeatureMat2 = [];
    
    % build feature matrix of this wordPair from all channel data
    [wordPairFeatureMat, ~] = buildChannelFeatureMat(channels, wordPairDir);
    
    % get random indices to downsize the comparisons made between pair events
    randIndices = randsample(size(wordPairFeatureMat, 1), 60);
    wordPairFeatureMat = wordPairFeatureMat(randIndices, :, :);
    
    %%- build events X features X time matrix with events being
    %%lined up to be compared
    for i=1:size(wordPairFeatureMat, 1),
        samePairFeatureMat1 = cat(1, samePairFeatureMat1, repmat(wordPairFeatureMat(i,:,:), length(i+1:size(wordPairFeatureMat, 1)), 1, 1));
        samePairFeatureMat2 = cat(1, samePairFeatureMat2, wordPairFeatureMat(i+1:end, :, :));
    end
    samePairFeatureMat1 = permute(samePairFeatureMat1, [1 3 2]);
    samePairFeatureMat2 = permute(samePairFeatureMat2, [1 3 2]);
end