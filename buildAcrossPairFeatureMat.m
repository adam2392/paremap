function [pairFeatureMat1, pairFeatureMat2] = buildAcrossPairFeatureMat(wordone, wordtwo, sessionBlockDir)
    pairFeatureMat1 = []; % initialize returned pair feature matrices
    pairFeatureMat2 = [];

    %%- 01: EXTRACT FILE DIRS FOR TWO VOCALIZED WORDS
    wordoneDir = fullfile(sessionBlockDir, wordone);
    filesone = dir(fullfile(wordoneDir));
    filesone = {filesone(3:end).name};
    wordtwoDir = fullfile(sessionBlockDir, wordtwo);
    filestwo = dir(fullfile(wordtwoDir));
    filestwo = {filestwo(3:end).name};

    %%- 02: EXTRACT CHANNEL FREQUENCY FEATURES 
    firstPairFeatureMat = [];
    secondPairFeatureMat = [];
    for iChan=1:length(filesone) % loop through all channels and get zscored power mat
        fileOnePath = fullfile(wordoneDir, filesone{iChan});
        dataOne = load(fileOnePath);
        dataOne = dataOne.data;
        fileTwoPath = fullfile(wordtwoDir, filestwo{iChan});
        dataTwo = load(fileTwoPath);
        dataTwo = dataTwo.data;

        if isempty(firstPairFeatureMat),
            firstPairFeatureMat = dataOne.powerMatZ;
            secondPairFeatureMat = dataTwo.powerMatZ;
        else
            firstPairFeatureMat = cat(2, firstPairFeatureMat, dataOne.powerMatZ);
            secondPairFeatureMat = cat(2, secondPairFeatureMat, dataTwo.powerMatZ);
        end
    end

    %%- 03: BUILD FEATURE MATRIX
    % build up the events X features X time matrix
    for i=1:size(firstPairFeatureMat, 1),
        pairFeatureMat1 = cat(1, pairFeatureMat1, repmat(firstPairFeatureMat(i,:,:), size(secondPairFeatureMat, 1), 1, 1));
        pairFeatureMat2 = cat(1, pairFeatureMat2, secondPairFeatureMat(:, :, :));
    end
end
    