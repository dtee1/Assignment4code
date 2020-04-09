function [xnew, count] = NR(t,x)
global G C b;
h = 0.3e-6;
bl = [0;0;0;0;25*sin(2*pi*10e3*(t));10];
bll = [0;0;0;0;25*sin(2*pi*10e3*(t+h));10];
A = G + (2/h).*C;
B = ((2/h).*C  - G)*x - [0;0;10e-15*(exp(x(3)/0.026)-1);0;0;0] + bl + bll;
phi = [1;1;1;1;1;1];
deltax = [1;1;1;1;1;1];
normdx = norm(deltax);
normphi = norm(phi);
count = 1;


while(normdx > 10e-12 && normphi > 10e-12 && count < 50)
    f = sparse([0;0;(10e-15*(exp(x(3)/0.026)-1));0;0;0]);
    phi = (A*x) + f - B;
    fprime = sparse([0 0 0 0 0 0;0 0 0 0 0 0;0 0 ((10e-15/0.026)*exp(x(3)/0.026)) 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0]);
    J = A + fprime;  
    deltax = J\-phi;
    x = x + deltax;
    normdx = norm(deltax);
    normphi = norm(phi);
    count = count + 1;
end
count  = count -1;
xnew = x;

end
