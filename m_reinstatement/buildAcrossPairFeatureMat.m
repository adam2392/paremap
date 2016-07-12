function [pairFeatureMat1, pairFeatureMat2] = buildAcrossPairFeatureMat(wordone, wordtwo, sessionFirstBlockDir, sessionSecondBlockDir)
    pairFeatureMat1 = []; % initialize returned pair feature matrices
    pairFeatureMat2 = [];
    
    if strfind(sessionFirstBlockDir, 'multitaper')
        multitaper = 1;
    else
        multitaper = 0;
    end
    
    
    %%- 01: EXTRACT FILE DIRS FOR TWO VOCALIZED WORDS
    wordoneDir = fullfile(sessionFirstBlockDir, wordone);
    filesone = dir(fullfile(wordoneDir));
    filesone = {filesone(3:end).name};
    wordtwoDir = fullfile(sessionSecondBlockDir, wordtwo);
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
        
        if multitaper == 1,
            addpath('./m_oldAnalysis_anovaANDsinglechannel/');

            % array of frequency bands
            freqBandAr(1).name    = 'delta';
            freqBandAr(1).rangeF  = [0 4];          %[2 4]
            freqBandAr(2).name    = 'theta';
            freqBandAr(2).rangeF  = [4 8];          %[4 8]
            freqBandAr(3).name    = 'alpha';
            freqBandAr(3).rangeF  = [8 16];         %[8 12]
            freqBandAr(4).name    = 'beta';
            freqBandAr(4).rangeF  = [16 32];        %[12 30]
            freqBandAr(5).name    = 'low gamma';
            freqBandAr(5).rangeF  = [32 80];        %[30 70]
            freqBandAr(6).name    = 'high gamma';
            freqBandAr(6).rangeF  = [80 160];       %[70 150]
            freqBandAr(7).name    = 'HFO';
            freqBandAr(7).rangeF  = [160 400];      %[150 400]

            % set the frequency bands to certain ranges for plotting
            for iFB=1:length(freqBandAr),
                freqBandAr(iFB).centerF = mean(freqBandAr(iFB).rangeF);
                %freqBandAr(iFB).label   = sprintf('%s-%.0fHz', freqBandAr(iFB).name(1:[ min( [length(freqBandAr(iFB).name), 6] )]), freqBandAr(iFB).centerF);
                freqBandAr(iFB).label   = sprintf('%s [%.0f-%.0f Hz]', freqBandAr(iFB).name, freqBandAr(iFB).rangeF);
            end
            freqBandYticks  = unique([freqBandAr(1:7).rangeF]);
            for iFB=1:length(freqBandYticks), freqBandYtickLabels{iFB} = sprintf('%.0f Hz', freqBandYticks(iFB)); end

            rangeFreqs = reshape([freqBandAr.rangeF], 2, 7)';
            freqs = squeeze(dataOne.freq);
            dataOne.powerMatZ = freqBinSpectrogram(dataOne.powerMatZ, rangeFreqs, freqs);
            dataTwo.powerMatZ = freqBinSpectrogram(dataTwo.powerMatZ, rangeFreqs, freqs);
        end
        
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
    