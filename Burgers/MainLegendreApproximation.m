global x Nburgers Nplot sols u xplot fappr expansioncoeff pols NQoI nuleft nuright N xburgers
f = @(x) burgers(0, x, Nburgers-1);
NxQ = N+1;
[x, w] = GaussLegendre(NxQ);
x = rescale(x, -1, 1, nuleft, nuright);
w = w./sum(w);
%sample solutions
sols = zeros(length(x), Nburgers);
for i = 1:length(x)
    [xburgers, sols(i, :)] = f(x(i));
end

u = sols(:, NQoI);

%% Compute Values of Polynomial Basis on the Quadrature Grid
pols(N+1).vals = [];
pols(N+1).pol = [];
for i = 1:N+1
    pols(i).pol = LegendrePol(i-1);
    pols(i).vals = polyval(pols(i).pol, rescale(x, nuleft, nuright, -1, 1));
end

expansioncoeff = zeros(N+1, 1);
for i = 1:N+1
    %Expansion Coefficients
    expansioncoeff(i) = DiscreteInner(w, u, pols(i).vals)/DiscreteInner(w, pols(i).vals, pols(i).vals);
end
%% Plot Approximation
xplot = linspace(nuleft, nuright, Nplot);
fappr = zeros(size(xplot));
for i = 1:N+1
    fappr = fappr + expansioncoeff(i)*polyval(pols(i).pol, rescale(xplot, nuleft, nuright, -1, 1));
end