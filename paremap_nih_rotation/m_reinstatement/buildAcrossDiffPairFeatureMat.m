function [pairFeatureMat1, pairFeatureMat2 ] = buildDiffPairFeatureMat(wordGroup, sessionFirstBlockDir, sessionSecondBlockDir)
    % loop through different words
    pairFeatureMat1 = [];
    pairFeatureMat2 = [];
    
    % loop through word groups (e.g. same, reverse, probe, different, target)
    for iPair=1:length(wordGroup),
        pairone = wordGroup{iPair}{1};
        pairtwo = wordGroup{iPair}{2};

        % create directories to the channel files for this word in this
        % block
        paironeDir = fullfile(sessionFirstBlockDir, pairone)
        paironeFiles = dir(fullfile(paironeDir));
        paironeFiles = {paironeFiles(3:end).name};
        
        pairtwoDir = fullfile(sessionSecondBlockDir, pairtwo)
        pairtwoFiles = dir(fullfile(pairtwoDir));
        pairtwoFiles = {pairtwoFiles(3:end).name};

        %%- Loop through channels of data
        firstPairFeatureMat = []; % eventsXfeaturesXtime
        secondPairFeatureMat = [];
        for iChan=1:length(paironeFiles) 
            paironeFilePath = fullfile(paironeDir, paironeFiles{iChan});
            pairtwoFilePath = fullfile(pairtwoDir, pairtwoFiles{iChan});

            dataPairOne = load(paironeFilePath);
            dataPairOne = dataPairOne.data;
            dataPairTwo = load(pairtwoFilePath);
            dataPairTwo = dataPairTwo.data;

            %%- concatenate all freq. vectors
            if isempty(firstPairFeatureMat)
                firstPairFeatureMat = dataPairOne.powerMatZ;
                secondPairFeatureMat = dataPairTwo.powerMatZ;
            else
                firstPairFeatureMat = cat(2, firstPairFeatureMat, dataPairOne.powerMatZ);
                secondPairFeatureMat = cat(2, secondPairFeatureMat, dataPairTwo.powerMatZ);
            end
        end

        %%- build events X features X time matrix
        for i=1:size(firstPairFeatureMat, 1)
            pairFeatureMat1 = cat(1, pairFeatureMat1, repmat(firstPairFeatureMat(i,:,:), size(secondPairFeatureMat, 1), 1, 1));
            pairFeatureMat2 = cat(1, pairFeatureMat2, secondPairFeatureMat(:, :, :));
        end
    end
end