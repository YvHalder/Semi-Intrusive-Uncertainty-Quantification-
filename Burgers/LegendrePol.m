function pol = LegendrePol(n)
if n == 0
    pol  = 1;
elseif n==1
    pol = [1; 0];
else
    pol1 = zeros(n+1, 1);
    pol2 = zeros(n+1, 1);
    pol2(end) = 1;
    pol1(end-1) = 1;
    for i = 2:n
        pol = ((2*i-1)/i)*circshift(pol1, -1) - ((i-1)/i)*pol2;
        pol2 = pol1;
        pol1 = pol;
    end
end