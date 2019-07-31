% Gauss-Seidel 

function [out,k] = GaussSeidel(x,y,dx,dy,RHS,p)

if(round(dx,3) ~= round(dy,3))
    disp('increments in x and y are not equivalent')
%     exit
end

M = size(x,1);
N = size(y,2);

tol = 1e-4;
k = 0;

pn = p;
pi = p;
% err = 1; % Can use if required
res = 1;
while res >tol
   
    k = k+1;
    
    for i = 2:M-1
        for j = 2:N-1     
        if i == M-1
            pn(i,j) = -(1/3).*(((dx^2).*RHS(i,j)) - (pi(i-1,j)+pi(i,j-1)+pi(i,j+1)));
        elseif j == 2
            pn(i,j) = -(1/3).*(((dx^2).*RHS(i,j)) - (pi(i-1,j)+pi(i,j+1)+pi(i+1,j)));
        else
        pn(i,j) = -0.25.*(((dx^2).*RHS(i,j)) - (pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
        end
        pi(i,j) = pn(i,j);
        end
    end
    pi(M,:) = pn(M-1,:);
    pi(:,1) = pn(:,2);
    pn(M,:) = pn(M-1,:);
    pn(:,1) = pn(:,2);
%     err = norm(pn-p,2);
%     p = pi;
    res = 0;
    for i = 2:M-1
        for j = 2:N-1
        res = res+ abs(((pn(i-1,j)+pn(i,j-1)+pn(i,j+1)+pn(i+1,j)-4*pn(i,j))./dx^2)-RHS(i,j));
        end
    end
    
end

out = pn;
end