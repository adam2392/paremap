function [diffPairFeatureMat1, diffPairFeatureMat2] = createDiffPairReinInput(channels, wordPairDir1, wordPairDir2)
    diffPairFeatureMat1 = [];
    diffPairFeatureMat2 = [];
    
    % build feature matrix of this wordPair from all channel data
    [wordPairFeatureMat1, wordPairFeatureMat2] = buildChannelFeatureMat(channels, wordPairDir1, wordPairDir2);
    
    % get random indices to downsize the comparisons made between pair events
    randIndices = randsample(min(size(wordPairFeatureMat1, 1), size(wordPairFeatureMat2, 1)), 60);
    wordPairFeatureMat1 = wordPairFeatureMat1(randIndices, :, :);
    wordPairFeatureMat2 = wordPairFeatureMat2(randIndices, :, :);
    
    %%- build events X features X time matrix with events being
    %%lined up to be compared
    for i=1:size(wordPairFeatureMat1, 1),
        diffPairFeatureMat1 = cat(1, diffPairFeatureMat1, repmat(wordPairFeatureMat1(i,:,:), size(wordPairFeatureMat2, 1), 1, 1));
        diffPairFeatureMat2 = cat(1, diffPairFeatureMat2, wordPairFeatureMat2(:, :, :));
    end
    diffPairFeatureMat1 = permute(diffPairFeatureMat1, [1 3 2]);
    diffPairFeatureMat2 = permute(diffPairFeatureMat2, [1 3 2]); 
end