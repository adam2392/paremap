close all
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/figures')
addpath('/Users/adam2392/Documents/MATLAB/SuPlot')

szindices = [61 74 81]; % for p04n007

% load saved figures
a = hgload('iiZone1_PY04N007.fig');
b = hgload('qAll74PY04N007.fig');
climit = get(gca, 'CLim');
c = hgload('cluster74All_PY04N007.fig');

xlimit = get(b, 'XLim')

% prepare subplots
figure
h(1) = subplot(2, 1, 1);
h(2) = subplot(2, 1, 2);
% h(3) = subplot(2, 2, 3);
% h(4) = subplot(2, 2, 4);

% paste figures on subplots
copyobj(allchild(get(b, 'CurrentAxes')), h(1));
copyobj(allchild(get(c, 'CurrentAxes')), h(2));
% copyobj(allchild(get(c, 'CurrentAxes')), h(4));

% modify their settings
% get(h(1))
% set(h(1), 'colormap', 'jet')
% set(h(1), 'YLabel', 'Channels');
% ax = axes;
% set(ax, 'YTick', [1 20 40 60 75])
% colormap('jet');
% c = colorbar;
% ylabel('channels');
% xlabel('time (seconds)');
% xlim([iitimes(j,1) iitimes(j,2)])
% title(['Interictal EVC W/O induced Seizures On zone: ' num2str(j)])
% set(gca, 'CLim', clim); 

% add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')

xlabel(h(1), 'time (seconds)');
xlabel(h(2), 'time (seconds)');
ylabel(h(1),'channels');
ylabel(h(2),'channels');
colormap(h(1), 'jet');
colorbar(h(1));
colorbar(h(2));
title(h(1), ['Pre -> Seizure -> Post ictal EVC and clustering for index 61'])
ylim(h(1), [0 75])
ylim(h(2), [0 10])
topxlim = get(h(1), 'XLim');
topxlim(1) = topxlim(1) + 20;
xlim(h(1), topxlim);
botxlim = get(h(2), 'XLim');
botxlim(1) = botxlim(1) + 20 ;
xlim(h(2), botxlim);