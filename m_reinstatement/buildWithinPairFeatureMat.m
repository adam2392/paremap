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
                freqs = squeeze(data.freq);
                data.powerMatZ = freqBinSpectrogram(data.powerMatZ, rangeFreqs, freqs);
            end
            
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
                
                %%- frequency bin the multitaper power for wordone and wordtwo
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
        
        %%- build up the events X features X time matrix
        for i=1:size(firstPairFeatureMat, 1),
            pairFeatureMat1 = cat(1, pairFeatureMat1, repmat(firstPairFeatureMat(i,:,:), size(secondPairFeatureMat, 1), 1, 1));
            pairFeatureMat2 = cat(1, pairFeatureMat2, secondPairFeatureMat(:, :, :));
        end
    end
end
    