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
           
            %%- check if this data comes from multitaper code
            if size(data.powerMatZ, 2) > 7, % need to condense along freq. domain first
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