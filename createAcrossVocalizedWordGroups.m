function wordPairs = createAcrossVocalizedWordGroups(targetWordsOne, targetWordsTwo)
    % all word pairs from the target words
    wordPairs = {};
    for iPair=1:length(targetWordsOne)
        firstWord = targetWordsOne{iPair}; % get current target word from block(i)

        for jPair=1:length(targetWordsTwo)
            secondWord = targetWordsTwo{jPair};
            
            % build up the word pairs list
            joinedWord = strjoin({firstWord, secondWord}, '_');
            wordPairs{end+1} = joinedWord;
        end
    end % loop through first target words
end