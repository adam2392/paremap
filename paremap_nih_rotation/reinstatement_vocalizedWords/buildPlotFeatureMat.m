function featureMat = buildFeatureMat(targetWord, sessionBlockDir)
    fileDir = fullfile(sessionBlockDir, targetWord);
    chanfiles = dir(fileDir);
    chanfiles = {chanfiles(3:end).name};
    
    %%- EXTRACT ALL CHANNEL DATA
    chanFile = fullfile(fileDir, chanfiles{1});
    data = load(chanFile);
    data = data.data;
    
    % initialize matrix to increase speed
    featureMat = zeros(size(data.powerMatZ, 1), ...
        size(data.powerMatZ, 2) * length(chanfiles), ...
        size(data.powerMatZ, 3));
    for iChan=1:length(chanfiles)
        chanFile = fullfile(fileDir, chanfiles{iChan});
        data = load(chanFile);
        data = data.data;
        
        %%- Concatenate all frequency vectors into feature vector
        featureMat(:, (iChan-1)*7+1:(iChan)*7, :) = data.powerMatZ;
    end
end