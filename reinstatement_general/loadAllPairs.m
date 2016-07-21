function [allWordPairs, allWordPairDirs] = loadAllPairs(dataDir)
    allWordPairs = dir(dataDir);
    allWordPairs = {allWordPairs(3:end).name};
    
    allWordPairDirs = fullfile(dataDir, allWordPairs);
end