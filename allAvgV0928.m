clc; clear vars; close all; tic;
addpath(genpath('D:\Papers and Data\Scitech 2020\Data\BOS\BOS_validation_supersonic_07232019\PIVmLibrary'))
MatLabSettings
load('ML.mat')

fb = 3;
ph = 80;
kgd = 0.238e-3;
%% Create list of images to analyze. 
numbers = ListPIVimages( img,{'AvgV'});
dt = 21171.8*1e-6;

    
    disp([num2str(1) ': ' strjoin(img(1).info) ' ...'])

    % Load image
    [x1,y1,U1,V1] = showimx(img(1).raw);
     vel1 = sqrt(U1.^2 + V1.^2);
    
     
   
    % Plot:
    figure(1)
    pcolor(x1,y1,vel1);
    shading flat
    colormap('jet(4096)')
    
%% Draw ROI (square) and then double click to proceed 
    % Note: there is a bug in the drawrectangle commands such that the size
    % in x and y is off by one. Redraw and try again in this case.
    
    h = drawrectangle(gca,'AspectRatio',1,'FixedAspectRatio',1);
    pos = customWait(h);
    
    
 %% Crop to selection
    xmin_pos = find(abs(x1(:,1)-pos(1,1))<1.5,1);
    xmax_pos = find(abs(x1(:,1)-pos(3,1))<1.5,1);
    ymax_pos = find(abs(y1(1,:)-pos(1,2))<1.5,1);
    ymin_pos = find(abs(y1(1,:)-pos(2,2))<1.5,1);
    
    x = x1(xmin_pos:xmax_pos,ymin_pos:ymax_pos)-x1(xmin_pos,1);
    y = y1(xmin_pos:xmax_pos,ymin_pos:ymax_pos)-y1(1,ymax_pos);
    vel = vel1(xmin_pos:xmax_pos,ymin_pos:ymax_pos);
    U = U1(xmin_pos:xmax_pos,ymin_pos:ymax_pos);
    V = V1(xmin_pos:xmax_pos,ymin_pos:ymax_pos);

    epx = (U*dt);
    epy = (V*dt);
     
    
    width = size(x,1);
    height = size(x,2);
    figure(2)
    pcolor(x,y,epx);
    shading interp
    colormap('jet(4096)')
    colorbar;
    
    
    %% Interpolate to modify grid spacing
    
    x = x./1000;
    y = y./1000;
    
    L = max(x(:,1));
    x = x./max(x(:,1));
    y = y./max(y(1,:));
    
    N = 100;

    x_interp = linspace(0,1,N);
    y_interp = fliplr(x_interp);
    
    xnew = ones(size(x_interp,2),size(x_interp,2));
    ynew = xnew;
    
    xnew = xnew.*x_interp';
    ynew = ynew.*y_interp;
    epx_new = interp2(x(:,1),y(1,:),epx,xnew(:,1),ynew(1,:));
    epy_new = interp2(x(:,1),y(1,:),epy,xnew(:,1),ynew(1,:));
    
    
    % Want to compare?

%     figure(5)
%     pcolor(x,y,epx); shading interp;
%     colorbar;
%     figure(6)
%     pcolor(xnew,ynew,epx_new); shading interp;
%     colorbar;
    
    %% Refractive index
    eta0 = 1.0029;
    rho_air = 1.204;
    zd = 394.335e-3; % Distance of background from schliere
    h = 72.54e-3; % Width of the schliere - taken from wedge dimensions
    M = 0.022; %scale factor * size of 1 pixel
 
    etax = eta0*L*epx_new/(M*zd*h*rho_air);
    etay = eta0*L*epy_new/(M*zd*h*rho_air);
    detax = zeros(N-2,N-2);
    detay = zeros(N-2,N-2);
    for i = 2:N-1
        for j = 2:N-1
            detax(i-1,j-1)=(etax(i+1,j)-etax(i-1,j))/2*(abs(xnew(2,1)-xnew(1,1)));
            detay(i-1,j-1)=(etay(i,j+1)-etay(i,j-1))/2*(abs(ynew(1,1)-ynew(1,2)));
        end
    end
    x2 = xnew(2:end-1,2:end-1);
    y2 = ynew(2:end-1,2:end-1);
    
    drho_x = etax/kgd;
    drho_y = etay/kgd;
    
    
%% Define RHS
   
  
    rhs = (detax+detay)./kgd;
    figure(3)
    pcolor(x2,y2,rhs);
    shading interp
    colormap('jet(4096)')
    colorbar;
    
    % for laplace equation set rhs to 0
    r1 = zeros(size(etax,1),size(etax,1));
    [rho,k,dx,dy] = Poisson_equation_2D(x2,y2,r1,r1);
    
    
%     [rho,k,dx,dy] = Poisson_equation_2D(x2,y2,etax/kgd,etay/kgd);
%     
    figure(4)
    pcolor(xnew,ynew,rho)
    
    shading interp
    colormap('jet(4096)')
    colorbar
%     set(gca,'clim',[1 1.7])
%     
%     figure(3)
%     pcolor(x(2:width-1,2:height-1),y(2:width-1,2:height-1),rho)
%     shading interp
%     colormap('jet')
%     
%     hold on;
%     quiver(x(2:width-1,2:height-1),y(2:width-1,2:height-1),detax./kgd,detay./kgd,'autoscale','on','autoscalefactor',5,'color','k')
%     xlim([-45 45])
%     ylim([-32 30])



