clear all;
clc;
close all;

disp('runnning')
brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'delta');
brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'theta');
brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'alpha');
brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'beta');
brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'low gamma');
brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'high gamma');
brainPlotSpectralPower('NIH034', 'multitaper', 'bipolar', 'HFO');

disp('nih037')
brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'delta');
brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'theta');
brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'alpha');
brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'beta');
brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'low gamma');
brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'high gamma');
brainPlotSpectralPower('NIH037', 'multitaper', 'bipolar', 'HFO');

brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'delta');
brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'theta');
brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'alpha');
brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'beta');
brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'low gamma');
brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'high gamma');
brainPlotSpectralPower('NIH039', 'multitaper', 'bipolar', 'HFO');