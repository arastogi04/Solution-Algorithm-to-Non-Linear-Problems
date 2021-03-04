clc;clear;
x = input("What is the value of x");
syms d;
d_value = zeros(40,1);
F1ext = linspace(0.25,10.0,40); %Load steps
i = 1;
for n = F1ext
    expr = (x*d/(10-d)) - 0.5*d*d;
    d_value(i) = eval(solve(expr==n,d,'Real',true));
    i = i+1;
end
plot(d_value,F1ext');
xlabel('d1');
ylabel('N1(d1)');
title("Exact Solution for x = "+x);