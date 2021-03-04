function G = LineSearchFunction (n,F2ext, d,del_d,s,x)

G = dot(del_d,Residual(n,F2ext,d+s*del_d,x));


end