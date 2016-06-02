function wordPairs = createWithinVocalizedWordGroups(targetWords)
    % all word pairs from the target words
    wordPairs = {};
    for i=1:length(targetWords)
        wordPairs{i} = {targetWords{i}, targetWords{i}};
    end
    allWordCombs = combnk(targetWords, 2)

    len = length(wordPairs);
    for i=1:length(allWordCombs)
        wordPairs{len+i} = {allWordCombs{i,1}, allWordCombs{i,2}};
    end
end