function plotEventSpectrogram(powerMat, binnedWaveT, freqBandAr, titleStr, subplotIndex, figSpect, NUMSUBPLOTS, TRIGGER)
    %% EXTRACT DATA FROM EVENTS
    %%%%%%%%%%%%%%%%%%%%%RECENTER THE TIMES WRT MSTIME %%%%%%%%%%%%%%%%
    powPlot = squeeze(mean(powerMat, 1));
    
    freqBandYtickLabels = {freqBandAr.label};
    %% BAR CHARTS
    figure(figSpect);
    subplot(NUMSUBPLOTS, 1,subplotIndex)
    
    hImg    = imagesc(binnedWaveT,[1:7],powPlot); 
    hold on;  colormap(jet);

    hCbar = colorbar('east');
    set(hCbar,'ycolor',[1 1 1]*.1, 'YAxisLocation', 'right')
            
    % set the heat map settings
    set(gca,'ytick',[1:7],'yticklabel',freqBandYtickLabels)
    set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
    title(sprintf('Plotting Events Power Matrix Z Scored for %s with %s events %s', ...
        TRIGGER, num2str(size(powerMat,1)), titleStr))

    % label the ranges for each band
    xlim=get(gca, 'xlim');
    rangeF = [1:7; 2:8]';
    for iFB=1:length(freqBandAr),
        hB = plot(xlim, [1 1]*rangeF(iFB,1)+0.5, '--', 'linewidth', 2, 'color', [1 1 1]*.7);
%         hB = plot(xlim, [1 1]*rangeF(iFB,2), '--', 'linewidth', 2, 'color', [1 1 1]*.7);
%         text(xlim(1)+range(xlim)*.02, mean(rangeF(iFB,:)) , freqBandAr(iFB).name, ...
%             'FontSize', 16, 'horizontalalignment','left');
    end
end