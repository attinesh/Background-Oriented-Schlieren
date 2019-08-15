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
res = zeros(M+1,M+2);
% Impose BC on pn

pn(1,:) = p(1,:);
% pn(M,:) = p(M,:);
done = false;

% Add ghost points on Neumann BC

pn = [pn; pn(end,:)]; % for i = M
pn = horzcat(pn(:,1), pn); % for j = 1
pn = horzcat(pn, pn(:,end)); % fpr j = M

pi = [pi; pi(end,:)]; % for i = M
pi = horzcat(pi(:,1), pi); % for j = 1
pi = horzcat(pi, pi(:,end)); % fpr j = M

dr_x = [dr_x; dr_x(end,:)]; % for i = M
dr_x = horzcat(dr_x(:,1), dr_x); % for j = 1
dr_x = horzcat(dr_x, dr_x(:,end)); % fpr j = M

dr_y = [dr_y; dr_y(end,:)]; % for i = M
dr_y = horzcat(dr_y(:,1), dr_y); % for j = 1
dr_y = horzcat(dr_y, dr_y(:,end)); % fpr j = M




while ~done
   
    k = k+1;
    
    for i = 2:M
        for j = 2:M+1                 
                pn(i,j) = -0.25.*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j)+dr_y(i,j+1)-dr_y(i,j-1))) - (pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
%                 pn(i,j) = 0.25.*((pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
                res(i,j) = pi(i,j) + 0.25.*(((dx/2).*(dr_x(i+1,j)-dr_x(i-1,j)+dr_y(i,j+1)-dr_y(i,j-1))) - (pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j)));
%                 res(i,j) = pi(i,j) - 0.25.*(pi(i-1,j)+pi(i,j-1)+pi(i,j+1)+pi(i+1,j));
        end
    end
    

    if norm(res(:),2) <tol
            done = true;
    end
    
pi = pn;
end
out = pn(1:M,2:M+1);
end

