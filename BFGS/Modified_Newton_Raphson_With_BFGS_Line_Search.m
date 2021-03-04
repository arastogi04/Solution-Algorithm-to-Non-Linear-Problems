clc;clear;
syms d1 d2;

x = input("What is the value of x");

N = [(x*d1/(10-d1))-0.5*d2*d2;(d2-d1)];

K = [diff(N(1),d1) diff(N(1),d2);diff(N(2),d1) diff(N(2),d2)];

F2ext = 0;

F1ext = linspace(0.25,10.0,40); %Load steps

r_tol = 1e-6; %Residual tolerance
itr_max = 15; %maximum newton-iterations
d = [0;0];
d_value = zeros(40,1);
count = 1;
V = zeros(20,2);
W = zeros(20,2);
for n = F1ext
    
    d_new = d;
    i = 1;
    K_init = eval(subs(K,[d1,d2],[d_new(1),d_new(2)]));
    G = 
    while i<itr_max
        if (i == 0)
            i = i+1;
            res = [n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)]));
            
            del_d = K_init\res;
            d = d+del_d;
            if (norm([n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)])))<r_tol)
                
                d_value(count) = d(1);
                
                count = count+1;
                
                break;
                
            end
        
        end
        del_res = eval(subs(N,[d1,d2],[d(1),d(2)])) - eval(subs(N,[d1,d2],[d(1)-del_d(1),d(2)-del_d(2)]));
        V(i,1:2) = del_d/(dot(del_d,del_res));
        alpha = sqrt(-1*dot(del_res,del_d)/(dot([n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)])),del_d))); 
        W(i,1:2) = -del_res + alpha*([n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)])));        
        
        for k = 1:i
            
        end
       
        res = [n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)]));
        del_d = eval(subs(K,[d1,d2],[d_new(1),d_new(2)]))\res;
        d = d+del_d;
        
        
        if (norm([n;F2ext] - eval(subs(N,[d1,d2],[d(1),d(2)])))<r_tol)
           
            d_value(count) = d(1);
            count = count+1;
            break;
            
        end
        
    end
    
end


plot(d_value,F1ext');
xlabel('d1');
ylabel('N1(d1)');
title("Modified Newton-Raphson for x = "+x);