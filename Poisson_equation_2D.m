
function [p,k,dx,dy]=Poisson_equation_2D(x,y,dr_x,dr_y)

% Solving the 2-D Poisson equation 

[Nx,Ny]=size(dr_x);


%Specifying parameters 

dx= abs(x(2,1)-x(1,1));            %Width of space step(x)
dy= abs(y(1,2)-y(1,1));            %Width of space step(y)
% x= x - min(x(:,1));                %Reset range from 0 to max
% y= y - min(y(1,:));                %Reset range from 0 to max
b=zeros(Nx,Ny);                    %Preallocating b


%%
% Initial Conditions
p=ones(Nx,Ny);                  %Preallocating p

%% 
% Rhs=fliplr(Rhs);   
% b=(Rhs);


%Boundary conditions (Note: Dirichlet should be normalized)

%% Top and Bottom
% Periodic boundary conditions

% p(:,1) = p(:,end);

% Dirichlet's conditions 
%     p(:,1)= 1; %top
%     p(:,end)= 1;% bottom
    
% Neumann's conditions
%     p(:,1)=p(:,2); %top
%     p(:,end)=p(:,end-1); % bottom

%% Left and Right

% Neumann's conditions   
%     p(1,:)=p(2,:); %Left
%     p(end,:)=p(end-1,:);%Right
    
% % Dirichlet's Conditions
    p(1,:)= 1; %Left
%     p(end,:)= 1.724;%Right

   
% Poisson equation solution (Based on specified method (SOR, Gauss-Seidel etc))
 [p,k] = Jacobi(x,y,dx,dy,dr_x,dr_y,p);
   
end			

 
