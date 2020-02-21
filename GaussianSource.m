function [A,x,y] = GaussianSource(m,sigma)
[xx,yy] = size(m);
x = linspace(-1,1,xx);
A = zeros(xx,yy);
g_scalar = 1/(sigma*sqrt(2*pi));
for i = 1:xx
    for j = 1:xx
g_exp = exp(-x(i).^2./(2.*sigma.^2));
y = g_scalar.*g_exp;
A(i,j)= y;
    end
end
end

