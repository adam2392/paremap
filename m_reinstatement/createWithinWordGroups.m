function [sameGroup, reverseGroup, diffGroup] = createWithinWordGroups(wordpairs)
    sameGroup = {};
    diffGroup = {};
    reverseGroup = {};
    
    %%- Create same word pairs group
    for i=1:length(wordpairs),
        sameGroup{i} = {wordpairs{i}, wordpairs{i}};
    end

    % find all combinations of words
    allWordCombs = combnk(wordpairs, 2);
    
    %%- Create reverse and different word pair group
    for i=1:length(allWordCombs),
        firstword = strsplit(allWordCombs{i,1}, '_');
        secondword = strsplit(allWordCombs{i,2}, '_');
        if strcmp(firstword(1), secondword(2)) && strcmp(firstword(2), secondword(1)),
            reverseGroup{end+1} = {allWordCombs{i,:}};
        else
            diffGroup{end+1} = {allWordCombs{i,:}};
        end
    end
end