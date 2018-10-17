% Generate numerical examples of Granger Causality

%% 01: 500 time points, 50 trials

xmodel = arima('Constant', 0, 'AR', {0.9, -0.5}, 'Variance', 0.1);
ymodel = arima('Constant', 0, 'AR', {0.8, -0.5}, 'Beta', [0.16, -0.2], 'Variance', 0.7);

ymodel

X = simulate(xmodel, 500);
Y = simulate(ymodel, 500);

CovXY = [1, 0.4; 0.4, 0.7];

figure(1)
plot(X)

figure(2)
plot(Y)

