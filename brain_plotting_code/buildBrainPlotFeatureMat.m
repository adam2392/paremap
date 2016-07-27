% return featureMat events X channels X time
function featureMat = buildFeatureMat(targetWord, sessionBlockDir, freqIndex)
    fileDir = fullfile(sessionBlockDir, targetWord);
    chanfiles = dir(fileDir);
    chanfiles = {chanfiles(3:end).name};
    
    %%- EXTRACT ALL CHANNEL DATA
    chanFile = fullfile(fileDir, chanfiles{1});
    data = load(chanFile);
    data = data.data;
    
    % initialize matrix to increase speed
    featureMat = zeros(size(data.powerMatZ, 1), ...
        length(chanfiles), ...
        size(data.powerMatZ, 3));
    for iChan=1:length(chanfiles)
        chanFile = fullfile(fileDir, chanfiles{iChan});
        data = load(chanFile);
        data = data.data;
        
        %%- Concatenate all frequency vectors into feature vector
        featureMat(:, iChan, :) = data.powerMatZ(:, freqIndex,:);
    end
end