% Returns:
% 1. allWordPairs: all the word pairs in the preprocessed directory
% 2. dataDir: the root directory to get access to each of the data in
% allWordPairs dir

function [allWordPairs, dataDir] = loadAllPairs(subj, TYPE_TRANSFORM, CUE_LOCK)
    dataDir = strcat('./condensed_data_', subj);
    dataDir = fullfile(dataDir, TYPE_TRANSFORM, CUE_LOCK)
    
    allWordPairs = dir(dataDir);
    allWordPairs = {allWordPairs(3:end).name};
end
