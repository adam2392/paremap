function [pairFeatureMat1, pairFeatureMat2] = buildWithinPairFeatureMat(wordone, wordtwo, sessionBlockDir)
    pairFeatureMat1 = []; % initialize returned pair feature matrices
    pairFeatureMat2 = [];
    
    if strfind(sessionBlockDir, 'multitaper')
        multitaper = 1;
    else
        multitaper = 0;
    end
    %%- if within blocks same pair feature mat
    if strcmp(wordone, wordtwo)
        disp('same')
        %%- 01: EXTRACT FILES FOR THE TWO WORD PAIRS
        wordDir = fullfile(sessionBlockDir, wordone);
        files = dir(fullfile(wordDir));
        files = {files(3:end).name};

        %%- 02: EXTRACT ALL CHANNEL DATA FOR EACH DIRECTORY
        wordPairFeatureMat = [];
        for iChan=1:length(files)
            filePath = fullfile(wordDir, files{iChan});
            data = load(filePath);
            data = data.data;

            %%- Concatenate all frequency vectors into feature vector
            if isempty(wordPairFeatureMat)
                wordPairFeatureMat = data.powerMatZ;
            else
                wordPairFeatureMat = cat(2, wordPairFeatureMat, data.powerMatZ);
            end 
        end

        %%- 03: BUILD FEATURE MATRIX
        % build up the events X features X time matrix
        for i=1:size(wordPairFeatureMat, 1),
            pairFeatureMat1 = cat(1, pairFeatureMat1, repmat(wordPairFeatureMat(i,:,:), length(i+1:size(wordPairFeatureMat, 1)), 1, 1));
            pairFeatureMat2 = cat(1, pairFeatureMat2, wordPairFeatureMat(i+1:end, :, :));
        end
    else
        disp('not same')
        
        %%- 01: EXTRACT FILES FOR TWO VOCALIZED WORDS
        wordoneDir = fullfile(sessionBlockDir, wordone);
        filesone = dir(fullfile(wordoneDir));
        filesone = {filesone(3:end).name};
        wordtwoDir = fullfile(sessionBlockDir, wordtwo);
        filestwo = dir(fullfile(wordtwoDir));
        filestwo = {filestwo(3:end).name};

        firstPairFeatureMat = [];
        secondPairFeatureMat = [];
        for iChan=1:length(filesone) % loop through all channels and get zscored power mat
            fileOnePath = fullfile(wordoneDir, filesone{iChan});
            dataOne = load(fileOnePath);
            dataOne = dataOne.data;
            fileTwoPath = fullfile(wordtwoDir, filestwo{iChan});
            dataTwo = load(fileTwoPath);
            dataTwo = dataTwo.data;
            
            %%- add the Z-scored power together as feature vectors
            if isempty(firstPairFeatureMat),
                firstPairFeatureMat = dataOne.powerMatZ;
                secondPairFeatureMat = dataTwo.powerMatZ;
            else
                firstPairFeatureMat = cat(2, firstPairFeatureMat, dataOne.powerMatZ);
                secondPairFeatureMat = cat(2, secondPairFeatureMat, dataTwo.powerMatZ);
            end
        end
        
        %%- build up the events X features X time matrix
        for i=1:size(firstPairFeatureMat, 1),
            pairFeatureMat1 = cat(1, pairFeatureMat1, repmat(firstPairFeatureMat(i,:,:), size(secondPairFeatureMat, 1), 1, 1));
            pairFeatureMat2 = cat(1, pairFeatureMat2, secondPairFeatureMat(:, :, :));
        end
    end
end
    