%% Semi-Intrusive Method for Simple ODE u_t = -alpha*u, with u(0)=1
%Exact solution: u(t) = exp(-alpha*t)
clear all
close all
clc

global x Nplot sols xplot fappr  pols al ar N f T fapprmin1 fapprplus1 fapprmin2 fapprplus2 expansioncoeff expansioncoeffmin1 expansioncoeffmin2 expansioncoeffplus1 expansioncoeffplus2
%We try to find an interpolation, which satisfies
f = @(a, t) exp(-a*t);

%Random Space: alpha lies in [al, ar]
al = 1; 
ar = 2;
%Uniform time mesh, which is used for determining the time derivatives
T = linspace(0, 1, 1000);

%% Legendre Approximation
N = 5; %Order of the Legendre Approximation (So we use N+1 quadrature points)
Nplot = 200; %The number of points where to evaluate the residue of the PDE (These points are uniformly spaced and we append the quadrature points as well)
MainLegendreApproximation; %Perform Legendre Approximation

%% Plotting
%Compute first Derivative (du/dt) and residue of the PDE in its original form (NOT including d/dalpha terms)
dtucentral4 = (fapprmin2 - 8*fapprmin1 + 8*fapprplus1 - fapprplus2)/(12*(T(end/2+1) - T(end/2))); %Fourth Order Central Differencing
dtucentral2 = (fapprplus1 - fapprmin1)/(2*(T(end/2+1) - T(end/2))); %Second Order Central Differencing
dtuforward = (fapprplus1 - fappr)/(T(end/2+1) - T(end/2)); %First order forward difference
dtubackward = (fappr - fapprmin1)/(T(end/2) - T(end/2-1)); %First order backward difference

residuecentral4 = xplot.*fappr + dtucentral4;
residuecentral2 = xplot.*fappr + dtucentral2;
residueforward = xplot.*fappr + dtuforward;
residuebackward = xplot.*fappr + dtubackward;

figure(1)
plot(xplot, abs(residuecentral4))
hold on
plot(xplot, abs(residuecentral2))
% plot(xplot, residueforward)
% plot(xplot, residuebackward)
plot(xplot, abs(fappr-f(xplot, T(end/2))))
plot(x, 1e-12.*ones(size(x)), 'o')
set(gca, 'YScale', 'log')
axis([al, ar, 1e-13, 1])
legend('Residue with Central Time Difference (fourth order)', 'Residue with Central Time Difference (second order)', 'Error Approximation: |fappr - fexc|')

%% Idea1: when, du/dt + alpha*u = Residue ==> du/dt = -alpha*u + R ==> du/dt = -alpha(u-R/alpha) ==> set u_new = u + coefficient(alpha)*Residue
% The coefficient is determined by sampling the solution at two other
% points for determining a linear trend for the term |fappr-fexc|/Residue
% The assumption from the linearity comes from the PDE itself

%DETERMINING COEFFICIENT(alpha) (linear)
samplelocations = [al; x(1:end-1) + diff(x)/2; ar]; %Sample Locations (end points of domain + in between quadrature nodes)
usamples = f(samplelocations, T(end/2)); %Exact Samples
%Approximate Samples
fsamp = zeros(size(samplelocations));
fsampmin1 = zeros(size(samplelocations));
fsampplus1 = zeros(size(samplelocations));
fsampmin2 = zeros(size(samplelocations));
fsampplus2 = zeros(size(samplelocations));
for i = 1:N+1
    fsamp = fsamp + expansioncoeff(i)*polyval(pols(i).pol, rescale(samplelocations, al, ar, -1, 1));
    fsampmin1 = fsampmin1 + expansioncoeffmin1(i)*polyval(pols(i).pol, rescale(samplelocations, al, ar, -1, 1));
    fsampplus1 = fsampplus1 + expansioncoeffplus1(i)*polyval(pols(i).pol, rescale(samplelocations, al, ar, -1, 1));
    fsampmin2 = fsampmin2 + expansioncoeffmin2(i)*polyval(pols(i).pol, rescale(samplelocations, al, ar, -1, 1));
    fsampplus2 = fsampplus2 + expansioncoeffplus2(i)*polyval(pols(i).pol, rescale(samplelocations, al, ar, -1, 1));
end
%Determine Residue at the sample locations
dtucentral4samp = (fsampmin2 - 8*fsampmin1 + 8*fsampplus1 - fsampplus2)/(12*(T(end/2+1) - T(end/2))); %Fourth Order Central Differencing
residuecentral4samp = samplelocations.*fsamp + dtucentral4samp;
coefficient = residuecentral4samp./(fsamp - usamples);
%Correct Approximation
fapprnew = fappr - 1./(coefficient(1) + (coefficient(end) - coefficient(1)).*rescale(xplot, al, ar, 0, 1)).*residuecentral4;
%Plot New Errors
figure(2)
plot(xplot, abs(fappr - f(xplot, T(end/2))), '-b')
hold on
plot(xplot, abs(fapprnew - f(xplot, T(end/2))), '-r')
set(gca, 'YScale', 'log')
axis([al, ar, 1e-13, 1])
legend('Old Error', 'New Error')



% %Compute Second Derivative (d^2u/(dalpha dt)) and residue formula 1
% %(including dalpha terms)
% dadtucentral = (derfapprplus1 - derfapprmin1)/(T(end/2+1) - T(end/2-1));
% dadtuforward = (derfapprplus1 - derfappr)/(T(end/2+1) - T(end/2));
% dadtubackward = (derfappr - derfapprmin1)/(T(end/2) - T(end/2-1));
% 
% residuecentral = derfappr + (1./fappr).*(fappr + dadtucentral);
% residueforward = derfappr + (1./fappr).*(fappr + dadtuforward);
% residuebackward = derfappr + (1./fappr).*(fappr + dadtubackward);

