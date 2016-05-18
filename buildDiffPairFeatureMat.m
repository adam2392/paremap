function [diffPairFeatureMat1, diffPairFeatureMat2 ] = buildDiffPairFeatureMat(diffWordGroup, sessionBlockDir)
    % loop through different words
    diffPairFeatureMat1 = [];
    diffPairFeatureMat2 = [];
    for iWord=1:length(diffWordGroup),
        wordone = diffWordGroup{iWord}{1};
        wordtwo = diffWordGroup{iWord}{2};

        %%- create directories to the channel files per word pair
        wordoneDir = fullfile(sessionBlockDir, wordone);
        filesone = dir(fullfile(wordoneDir));
        filesone = {filesone(3:end).name};
        wordtwoDir = fullfile(sessionBlockDir, wordtwo);
        filestwo = dir(fullfile(wordtwoDir));
        filestwo = {filestwo(3:end).name};

        %%- loop through channels of data
        firstPairFeatureMat = [];
        secondPairFeatureMat = [];
        for iChan=1:length(filesone) % loop through and open up all channels
            % Load in the data struct for each word pair per channel
            fileOnePath = fullfile(wordoneDir, filesone{iChan});
            dataOne = load(fileOnePath);
            dataOne = dataOne.data;
            fileTwoPath = fullfile(wordtwoDir, filestwo{iChan});
            dataTwo = load(fileTwoPath);
            dataTwo = dataTwo.data;

            dataOne.powerMatZ = dataOne.powerMatZ(:,:,1:size(dataOne.powerMatZ,3)-1);
            dataTwo.powerMatZ = dataTwo.powerMatZ(:,:,1:size(dataTwo.powerMatZ,3)-1);
            
            
            % concatenate all the freq. vectors that are already 500 ms
            % windowed and 100 ms overlap
            if isempty(firstPairFeatureMat),
                firstPairFeatureMat = dataOne.powerMatZ;
                secondPairFeatureMat = dataTwo.powerMatZ;
            else
                firstPairFeatureMat = cat(2, firstPairFeatureMat, dataOne.powerMatZ);
                secondPairFeatureMat = cat(2, secondPairFeatureMat, dataTwo.powerMatZ);
            end
        end
%             disp('different pair feature mats')
%             size(firstPairFeatureMat)
%             size(secondPairFeatureMat)
        %%- build events X features X time matrix with events being
        %%lined up to be compared
        for i=1:size(firstPairFeatureMat, 1),
            diffPairFeatureMat1 = cat(1, diffPairFeatureMat1, repmat(firstPairFeatureMat(i,:,:), size(secondPairFeatureMat, 1), 1, 1));
            diffPairFeatureMat2 = cat(1, diffPairFeatureMat2, secondPairFeatureMat(:, :, :));
        end
    end % end of loop through different words
end