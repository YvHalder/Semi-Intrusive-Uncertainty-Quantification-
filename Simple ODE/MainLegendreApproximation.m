global x Nplot sols xplot fappr  pols al ar N f T fapprmin1 fapprplus1 fapprmin2 fapprplus2 expansioncoeff expansioncoeffmin1 expansioncoeffmin2 expansioncoeffplus1 expansioncoeffplus2
NxQ = N+1;
[x, w] = GaussLegendre(NxQ);
x = rescale(x, -1, 1, al, ar);
w = w./sum(w);
%sample solutions
sols = f(x, T);

%Solution at different timelevels, used for determining the time derivatives
umin2 = sols(:, end/2-2);
umin1 = sols(:, end/2-1);
u = sols(:, end/2);
uplus1 = sols(:, end/2+1);
uplus2 = sols(:, end/2+2);

%% Compute Values of Polynomial Basis on the Quadrature Grid
pols(N+1).vals = [];
pols(N+1).pol = [];
for i = 1:N+1
    pols(i).pol = LegendrePol(i-1);
    pols(i).vals = polyval(pols(i).pol, rescale(x, al, ar, -1, 1));
end

expansioncoeff = zeros(N+1, 1);
expansioncoeffmin1 = zeros(N+1, 1);
expansioncoeffplus1 = zeros(N+1, 1);
expansioncoeffmin2 = zeros(N+1, 1);
expansioncoeffplus2 = zeros(N+1, 1);
for i = 1:N+1
    %Expansion Coefficients
    expansioncoeff(i) = DiscreteInner(w, u, pols(i).vals)/DiscreteInner(w, pols(i).vals, pols(i).vals);
    expansioncoeffmin1(i) = DiscreteInner(w, umin1, pols(i).vals)/DiscreteInner(w, pols(i).vals, pols(i).vals);
    expansioncoeffplus1(i) = DiscreteInner(w, uplus1, pols(i).vals)/DiscreteInner(w, pols(i).vals, pols(i).vals);
    expansioncoeffmin2(i) = DiscreteInner(w, umin2, pols(i).vals)/DiscreteInner(w, pols(i).vals, pols(i).vals);
    expansioncoeffplus2(i) = DiscreteInner(w, uplus2, pols(i).vals)/DiscreteInner(w, pols(i).vals, pols(i).vals);
end
%% Determine Approximation
xplot = linspace(al, ar, Nplot); %grid for determining the approximation/residual/exact solution
xplot = sort([xplot, x']); %append quadrature nodes

fappr = zeros(size(xplot));
fapprmin1 = zeros(size(xplot));
fapprplus1 = zeros(size(xplot));
fapprmin2 = zeros(size(xplot));
fapprplus2 = zeros(size(xplot));
for i = 1:N+1
    fappr = fappr + expansioncoeff(i)*polyval(pols(i).pol, rescale(xplot, al, ar, -1, 1));
    fapprmin1 = fapprmin1 + expansioncoeffmin1(i)*polyval(pols(i).pol, rescale(xplot, al, ar, -1, 1));
    fapprplus1 = fapprplus1 + expansioncoeffplus1(i)*polyval(pols(i).pol, rescale(xplot, al, ar, -1, 1));
    fapprmin2 = fapprmin2 + expansioncoeffmin2(i)*polyval(pols(i).pol, rescale(xplot, al, ar, -1, 1));
    fapprplus2 = fapprplus2 + expansioncoeffplus2(i)*polyval(pols(i).pol, rescale(xplot, al, ar, -1, 1));
end