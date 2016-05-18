function [samePairFeatureMat1, samePairFeatureMat2] = buildSamePairFeatureMat(sameWordGroup, sessionBlockDir)
    % loop through same word group
    samePairFeatureMat1 = [];
    samePairFeatureMat2 = [];
    
    % loop through every word in the sameWordGroup
    for iWord=1:length(sameWordGroup),
        word = sameWordGroup{iWord}{1};

        wordDir = fullfile(sessionBlockDir, word);
        files = dir(fullfile(wordDir));
        files = {files(3:end).name};

        wordPairFeatureMat = [];
        for iChan=1:length(files) % loop through and open up all channels
            filePath = fullfile(wordDir, files{iChan});
            data = load(filePath);
            data = data.data;
            data.powerMatZ = data.powerMatZ(:,:,1:size(data.powerMatZ,3)-1);
            
            % concatenate all the freq. vectors that are already 500 ms
            % windowed and 100 ms overlap
            if isempty(wordPairFeatureMat),
                wordPairFeatureMat = data.powerMatZ;
            else
                wordPairFeatureMat = cat(2, wordPairFeatureMat, data.powerMatZ);
            end
        end
%             disp('size of wordpairfeaturemat')
%             size(wordPairFeatureMat)
        %%- build events X features X time matrix with events being
        %%lined up to be compared
        for i=1:size(wordPairFeatureMat, 1),
            samePairFeatureMat1 = cat(1, samePairFeatureMat1, repmat(wordPairFeatureMat(i,:,:), length(i+1:size(wordPairFeatureMat, 1)), 1, 1));
            samePairFeatureMat2 = cat(1, samePairFeatureMat2, wordPairFeatureMat(i+1:end, :, :));
        end
    end
end