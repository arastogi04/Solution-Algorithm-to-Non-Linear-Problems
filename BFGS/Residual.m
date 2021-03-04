function res = Residual(n,F2ext,d,x)

syms d1 d2 s;
N = [(x*d1/(10-d1))-0.5*d2*d2;(d2-d1)];

res = [n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)]));

end