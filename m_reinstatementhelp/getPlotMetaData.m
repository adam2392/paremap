function [ticks, labels, timeZero] = getPlotMetaData(dataFile)
    exData = load(dataFile);
    data = exData.data;
    timeTicks = data.waveT(:,2);
    ticks = 1:5:length(timeTicks);
    labels = timeTicks(ticks);
    timeZero = data.timeZero;
end