subj = 'NIH034';

workingDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation';
subjDir = fullfile(workingDir, strcat('condensed_data_', subj), 'multitaper_global', 'vocalization');

sessionDirs = dir(subjDir);
sessionDirs = sessionDirs(3:end)
%%- loop over subjDir to get sessions
for iDir=1:length(sessionDirs)
    sessionDir = fullfile(subjDir, sessionDirs(iDir).name);
    blockDirs = dir(sessionDir);
    blockDirs = blockDirs(3:end);
    
    %%- loop over block directories
    for iBlock=1:length(blockDirs)
        blockDir = fullfile(sessionDir, blockDirs(iBlock).name);
        pairDirs = dir(blockDir);
        pairDirs = pairDirs(3:end);
        
        %%- loop over word pair directories
        for iPair=1:length(pairDirs)
            pairDir = fullfile(blockDir, pairDirs(iPair).name)
            
            %%- get all the channels and 
            chanMats = dir(pairDir);
            chanMats = chanMats(3:end);
            
            %%- loop over every channel mat file
            for iChan=1:length(chanMats)
                chanName = chanMats(iChan).name;
                chanMat = fullfile(pairDir, chanMats(iChan).name);
                
                % load in the data and condense it's frequency
                data = load(chanMat);
                data = data.data;
                
                %%- check if this data comes from multitaper code; if so,
                %%condense...
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
                else
                    disp('Didnt change');
                    disp(chanMat)
                end
                save(chanMat, 'data');
            end
        end
    end
end

