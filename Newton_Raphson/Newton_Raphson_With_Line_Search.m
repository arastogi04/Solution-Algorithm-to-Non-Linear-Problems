%Case (ii) - Newton-Raphson Method with line search

clc;clear;
syms d1 d2 s;

x = input("What is the value of x");

N = [(x*d1/(10-d1))-0.5*d2*d2;(d2-d1)]; %F_internal

K = [diff(N(1),d1) diff(N(1),d2);diff(N(2),d1) diff(N(2),d2)]; %Consistent Tangent

F2ext = 0;
iterations = zeros(40,1);
F1ext = linspace(0.25,10.0,40); %Load steps

r_tol = 1e-6; %Residual tolerance
itr_max = 15; %maximum newton-iterations
d = [0;0];
d_value = zeros(40,1);
count = 1;

for n = F1ext   %Load step loop
    
    i = 0;
    while i<itr_max %Newton iterations loop
        
        i = i+1;
        res = [n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)])); %Residual at iteration i
        
        del_d = eval(subs(K,[d1,d2],[d(1),d(2)]))\res; %Solving the linear equation for del_d at i
        
        d_temp = d+s*del_d;
        
        expr = dot(del_d,[n;F2ext] - eval(subs(N,[d1,d2],[d_temp(1),d_temp(2)])));
        s_value = eval(solve(expr==0,s,'Real',true));
        d = d+s_value*del_d; %Update the solution vector for the iteration i+1 (stop if solution converged)
        
        if (norm([n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)])))<r_tol) %Check the convergence
            
            iterations(count) = i;
            d_value(count) = d(1);
            count = count+1;
            break;
        end
        
    end
    
end

plot(F1ext',iterations);
xlabel('Load Steps');
ylabel('Iterations');
title("Newton-Raphson with line search for x = "+x);

plot(d_value,F1ext');
xlabel('d1');
ylabel('N1(d1)');
title("Numerical Solution (Newton-Raphson with Line Search) for x = "+x);