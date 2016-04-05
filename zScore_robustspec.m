%%% Script to Z score all robust spect output
condensedDataDir = '/Users/adam2392/Documents/MATLAB/Johns Hopkins/NINDS_Rotation/condensed_data/robust_spec/';

ext = '*.mat';
files = dir(fullfile(condensedDataDir, ext));
files = {files.name};
file = strcat(condensedDataDir, files{1});
data = load(file);
data = data.data;

for i=1:length(files)
    file = fullfile(condensedDataDir, files{i});
    data = load(file);
    data = data.data;
    
    powerMatZ = data.powerMatZ;
    % log the power matrix first
    powerMatZ = 10*log10(powerMatZ);
    
    % initialize new power matrix to be z-scored
    newpowerMatZ = zeros(size(powerMatZ));
    
    % for each file, z-score each channel and frequency
    tic;
    fprintf(' [%.1f sec] --> z-score', toc);  tic;
    for iF = 1:length(size(powerMatZ,2)),
        allVal = squeeze(powerMatZ(:,iF,:));
        allVal = squeeze(allVal(:)); %allVal for particular chan and freq
        mu = mean(allVal); stdev = std(allVal);

        % create the power matrix
        newpowerMatZ(:,iF,:) = (powerMatZ(:,iF,:)-mu)/stdev;

        if sum(isnan(newpowerMatZ(:,iF,:)))>0
            keyboard;
        end
    end
    fprintf(' [%.1f sec]', toc); tic;
    
    %%- Save the newpowermatrix into the data struct
    
end