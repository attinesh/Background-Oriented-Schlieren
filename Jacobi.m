% Jacobi 

function [out,k] = Jacobi(x,y,dx,dy,dr_x,dr_y,p)

if(round(dx,3) ~= round(dy,3))
    disp('increments in x and y are not equivalent')
%     exit
end

M = size(x,1);
N = size(y,2);

tol = 1e-5;
k = 0;

pn = p;
pi = p;
% err = 1; % Can use if required
res = 1;
while res >tol
   
    k = k+1;
    
    for i = 2:M-1
        for j = 1:N     
            if j == N 
                pn(i,j) = -(1/4).*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j))) - (pi(i-1,j)+pi(i+1,j)+ 2*pi(i,j-1)));
            elseif j == 1
                pn(i,j) = -(1/4).*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j))) - (pi(i-1,j)+2*pi(i,j+1)+pi(i+1,j)));
 
            else
                pn(i,j) = -0.25.*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j)+dr_y(i,j+1)-dr_y(i,j-1))) - (pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
            end
        end
    end


    res = 0;
    for i = 2:M-1
        for j = 2:N-1
        res = res+ abs(((pn(i-1,j)+pn(i,j-1)+pn(i,j+1)+pn(i+1,j)-4*pn(i,j))./dx^2)-((dr_x(i+1,j)-dr_x(i-1,j)+dr_y(i,j+1)-dr_y(i,j-1))./(2*dx)));
        end
    end
pi = pn;   
end

out = pn;
end