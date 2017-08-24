global x Nburgers Nplot sols u xplot fappr expansioncoeff pols NQoI nuleft nuright N xburgers x2 u2 sols2
f = @(x) burgers(0, x, Nburgers-1);
x2 = linspace(nuleft, nuright, N+1);

% Sample solutions
sols2 = zeros(length(x2), Nburgers);
for i = 1:length(x)
    [xburgers, sols2(i, :)] = f(x2(i));
end

u2 = sols2(:, NQoI);

