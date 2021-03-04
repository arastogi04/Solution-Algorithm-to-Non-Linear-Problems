clc;clear;
syms d1 d2 s;

x = input("What is the value of x");

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
    
    i = 1; 
     
    K_init = Tangent(d,x);    
    
    del_d = K_init\Residual(n,F2ext,d,x);
    
    %G_expr = dot(del_d,Residual(n,F2ext,d+s*del_d,x));
    
    %s_value = eval(solve(G_expr==0,s,'Real',true));
    s_value = 1;
    [v,w] = BFGS_V_W(n,F2ext,d,del_d,s_value,x);
    
    V(1,1:2) = v';
    W(1,1:2) = w';
    
    d = d+s_value*del_d;
    
    R_bar = Residual(n,F2ext,d,x)';
    
    i = i+1;
    
    while i<itr_max                
        
        
        for k = 1:i-1
            
            R_bar = R_bar*(eye(2,2) + W(k,1:2)' * V(k,1:2));
            
        end
        
        del_d_bar = (K_init\R_bar')';
        
        for k = 1:i-1
            
            del_d_bar = del_d_bar * (eye(2,2) + V(k,1:2)'*W(k,1:2));
            
        end
        
        del_d = del_d_bar';
        
        %G_expr = dot(del_d,Residual(n,F2ext,d+s*del_d,x));
    
        %s_value = eval(solve(G_expr==0,s,'Real',true));
        [v,w] = BFGS_V_W(n,F2ext,d,del_d,s_value,x);
        V(i,1:2) = v';
        W(i,1:2) = w';
        i = i+1;
        d = d+s_value*del_d;
        disp(i);
        disp(norm(Residual(n,F2ext,d,x)));
        if (norm(Residual(n,F2ext,d,x))<r_tol)
           
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