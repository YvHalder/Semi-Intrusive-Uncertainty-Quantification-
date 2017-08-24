%% Semi-Intrusive Burgers Equation
warning('off')
clc
clear all
global x Nburgers Nplot sols u xplot fappr expansioncoeff pols NQoI nuleft nuright N xburgers x2 u2 sols2

%% Legendre Approximation
Nburgers = 10000;
NQoI = floor(1*Nburgers/4);
nuleft = 0.1;
nuright = 1;
N = 1000;
Nplot = 200;

MainLegendreApproximation;
MainUniformSampling;

%% Determine Derivatives
firstder = (sols2(:, NQoI+1) - sols2(:, NQoI-1))./(2*(xburgers(10)-xburgers(9)));
secder = (sols2(:, NQoI-1) - 2*sols2(:, NQoI) + sols2(:, NQoI+1))./((xburgers(10)-xburgers(9)).^2);

%Term 1
term1 = secder;

%Term2
dernu = (secder(2:end) - secder(1:end-1))./(x2(2:end) - x2(1:end-1))';
dernu(end+1) = dernu(end);
term2 = x2'.*dernu;

%Term3
dernu = (firstder(2:end) - firstder(1:end-1))./(x2(2:end) - x2(1:end-1))';
dernu(end+1) = dernu(end);
term3 = sols2(:, NQoI).*dernu;

rhs = (term1+term2-term3)./firstder;

plot(x(1:end-1), diff(u)./diff(x))
hold on
plot(x2(1:end-1), diff(u2)./diff(x2))
plot(x2, rhs, '-s')