function [sameGroup, reverseGroup, probeGroup, targetGroup, diffGroup] = createWithinWordGroups(firstwordpairs, secondwordpairs)
    DEBUG_ON = 0;

    sameGroup = {};
    reverseGroup = {};
    probeGroup = {};
    targetGroup = {};
    diffGroup = {};

    %%- Loop through first wordpairs and compare with block(i+1)
    for iPair=1:length(firstwordpairs),
        firstWordPair = firstwordpairs{iPair}; % get the current word pair
        
        %%- Get the indices for different groupings in the next blocks'
        %%word pairs (secondwordpairs)
        sameWordIndices = findSameIndices(firstWordPair, secondwordpairs);
        reverseWordIndices = findReverse(firstWordPair, secondwordpairs);
        diffWordIndices = findDifferent(firstWordPair, secondwordpairs);
        probeWordIndices = findProbe(firstWordPair, secondwordpairs);
        targetWordIndices = findTarget(firstWordPair, secondwordpairs);
        
        
        if DEBUG_ON
            disp('Checking same word find.')
            firstWordPair
            secondwordpairs{sameWordIndices}
            disp('checking reverse word find.')
            firstWordPair
            secondwordpairs{reverseWordIndices}
            disp('checking diff word find.')
            firstWordPair
            secondwordpairs{diffWordIndices}
            firstWordPair
            disp('probe')
            secondwordpairs{probeWordIndices}
            disp('target')
            secondwordpairs{targetWordIndices}
        end
        
        %%- Add wordPairs to their corresponding group sets
        if (~isempty(sameWordIndices)) % check if we found any in secondwordgroup
            sameGroup{end+1} = {firstWordPair, secondwordpairs{sameWordIndices}};
        end
        if (~isempty(reverseWordIndices))
            reverseGroup{end+1} = {firstWordPair, secondwordpairs{reverseWordIndices}};
        end
        if (~isempty(diffWordIndices))
            % diffWordIndices could be a list instead of just 1 pairing
            if length(diffWordIndices) == 1
                diffGroup{end+1} = {firstWordPair, secondwordpairs{diffWordIndices}};
            else
                for diffIndex=1:length(diffWordIndices)
                    diffGroup{end+1} = {firstWordPair, secondwordpairs{diffWordIndices(diffIndex)}};
                end
            end
        end
        if (~isempty(probeWordIndices))
            % probeWordIndices could be a list instead of just 1 pairing
            if length(probeWordIndices) == 1
                probeGroup{end+1} = {firstWordPair, secondwordpairs{probeWordIndices}};
            else
                for probeIndex=1:length(probeWordIndices)
                    probeGroup{end+1} = {firstWordPair, secondwordpairs{probeWordIndices(probeIndex)}};
                end
            end
        end
        if (~isempty(targetWordIndices))
            % diffWordIndices could be a list instead of just 1 pairing
            if length(targetWordIndices) == 1
                targetGroup{end+1} = {firstWordPair, secondwordpairs{targetWordIndices}};
            else
                for targetIndex=1:length(targetWordIndices)
                    targetGroup{end+1} = {firstWordPair, secondwordpairs{targetWordIndices(targetIndex)}};
                end
            end
        end
    end % end of loop through word groups
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ HELPER FUNCTIONS ---------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sameWordIndices = findSameIndices(wordPair, otherWordPairGroup)
    words = strsplit(wordPair, '_');
    sameWordIndices = find(ismember(otherWordPairGroup, wordPair));
end
      
% find indices that show reversal of wordPair in wordPairComb
function reverseWordIndices = findReverse(wordPair, otherWordPairGroup)
    reverseWordPair = fliplr(strsplit(wordPair, '_')); % split word and reverse
    reverseWordPair = strjoin(reverseWordPair, '_');   % join the reversed word
    reverseWordIndices = find(ismember(otherWordPairGroup, reverseWordPair)); % find in cell array
end

% find indices that show different of wordPair in wordPairComb
function diffWordIndices = findDifferent(wordPair, otherWordPairGroup)
    words = strsplit(wordPair, '_'); % split word and reverse
    diffWordIndices = [];

    % loop through all wordpairs in other group
    for iPair=1:length(otherWordPairGroup),
        otherWordPair = strsplit(otherWordPairGroup{iPair}, '_');
        
        % append index, if the word is different 
        if ~ismember(words, otherWordPair)
            diffWordIndices = [diffWordIndices; iPair];
        end 
    end
end

function probeWordIndices = findProbe(wordPair, otherWordPairGroup)
    words = strsplit(wordPair, '_'); % split word and reverse    
    probeWordIndices = [];

    for iPair=1:length(otherWordPairGroup),
        otherWordPair = strsplit(otherWordPairGroup{iPair}, '_');

        % check that is is not the samepair and the probe word overlaps
        if isempty(find(ismember(otherWordPairGroup, wordPair))) && ...
                strcmp(words{1}, otherWordPair{1})

            probeWordIndices = [probeWordIndices; iPair];
        end 
    end
end

function targetWordIndices = findTarget(wordPair, otherWordPairGroup)
    words = strsplit(wordPair, '_'); % split word and reverse    
    targetWordIndices = [];

    for iPair=1:length(otherWordPairGroup),
        otherWordPair = strsplit(otherWordPairGroup{iPair}, '_');

        % check that is is not the samepair and the probe word overlaps
        if isempty(find(ismember(otherWordPairGroup, wordPair))) && ...
                strcmp(words{2}, otherWordPair{2})

            targetWordIndices = [targetWordIndices; iPair];
        end 
    end
end

