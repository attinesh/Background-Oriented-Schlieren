% Jacobi 

function [out,k] = Jacobi(x,y,dx,dy,dr_x,dr_y,p)

if(round(dx,3) ~= round(dy,3))
    disp('increments in x and y are not equivalent')
%     exit
end

M = size(p,1);


tol = 1e-5;
k = 0;

pn = zeros(M,M);
pi = p;
res = pn;
% Impose BC on pn

pn(1,:) = p(1,:);
pn(M,:) = p(M,:);
done = false;

while ~done
   
    k = k+1;
    
    for i = 2:M-1
        for j = 1:M     
            if j == M 
                pn(i,j) = -(1/4).*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j))) - (pi(i-1,j)+pi(i+1,j)+ 2*pi(i,j-1)));
%                  pn(i,j) = (1/4).*((pi(i-1,j)+pi(i+1,j)+ 2*pi(i,j-1))); 
                 res(i,j) = pi(i,j) + 0.25.*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j))) - (pi(i-1,j)+pi(i+1,j)+ 2*pi(i,j-1)));
            elseif j == 1
                pn(i,j) = -(1/4).*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j))) - (pi(i-1,j)+pi(i+1,j)+ 2*pi(i,j+1)));
%                 pn(i,j) = (1/4).*((pi(i-1,j)+pi(i+1,j)+ 2*pi(i,j+1)));
                res(i,j) = pi(i,j) + 0.25.*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j))) - (pi(i-1,j)+pi(i+1,j)+ 2*pi(i,j+1)));
 
            else
                pn(i,j) = -0.25.*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j)+dr_y(i,j+1)-dr_y(i,j-1))) - (pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
%                 pn(i,j) = 0.25.*((pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
                res(i,j) = pi(i,j) + 0.25.*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j)+dr_y(i,j+1)-dr_y(i,j-1))) - (pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
            end
        end
    end

    if norm(res(:),2) <tol
            done = true;
    end
    
pi = pn;   
end

out = pn;
end