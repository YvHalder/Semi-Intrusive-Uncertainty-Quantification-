function y = DiscreteInner(w, f, g)
% Computer the discrete innerproduct at the nodes x with corresponding
% weights w, of the function f and g.

y = sum(w.*f.*g);