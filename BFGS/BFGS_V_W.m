function [v,w] = BFGS_V_W (n,F2ext,d,del_d,s,x)

v = del_d/(LineSearchFunction (n,F2ext, d,del_d,s,x) - LineSearchFunction (n,F2ext, d,del_d,0,x));

alpha = sqrt((-s*(LineSearchFunction (n,F2ext, d,del_d,s,x) - LineSearchFunction (n,F2ext, d,del_d,0,x)))/LineSearchFunction (n,F2ext, d,del_d,0,x));

w = -1*(Residual(n,F2ext,d+s*del_d,x) - Residual(n,F2ext,d,x)) +  alpha*Residual(n,F2ext,d,x);

end