function K = Tangent(d,x)

syms d1 d2;
N = [(x*d1/(10-d1))-0.5*d2*d2;(d2-d1)];

K = [diff(N(1),d1) diff(N(1),d2);diff(N(2),d1) diff(N(2),d2)];
K = eval(subs(K,[d1,d2],[d(1),d(2)]));
end