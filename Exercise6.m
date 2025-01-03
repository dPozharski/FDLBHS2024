%%%%%%
%Exercise #6 Lattice Boltzman Methods
%Denis Pozharskiy
%23/12/2024
%%%%%%
clear;
close all;
setGroot();
colormap(redblue);

%Initialization and Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7  8  9
%  \ | /
%   \|/
% 6--1--2
%   /|\
%  / | \
% 5  4  3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NL = 9;
cx = [0 0 1 1 1 0 -1 -1 -1];
cy = [0 1 1 0 -1 -1 -1 0 1];
w = [4/9 1/9 1/36 1/9 1/36 1/9 1/36 1/9 1/36];
opp = [1,6,7,8,9,2,3,4,5];%Index of opposite population
dt =  1;
dr = 1;
Ds = [12,16,20,24,28,32,36];
Re = 20;
rho0 = 1;
T0 = 1/3;
cs = sqrt(T0);
U0 = 0.1;
maxT = 20000;
recirculationDist = zeros(3,1);

for iD = 1:length(Ds)
    D = Ds(iD);
    y = 0:22*D;
    x = 0:4.1*D;
    [X,Y] = ndgrid(x,y);
    Ux = zeros(size(X));
    Uy = zeros(size(X));

    rho = ones(size(X));
    nu = U0*D/Re;
    beta = (2*nu/T0 + 1)^-1;
    center = [2*D,2*D]; %Center of cyllindrical obstacle
    H = 4.1*D;
    dist = sqrt((X-center(1)).^2 + (Y-center(2)).^2);
    boundary = find(~(dist>D/2)); %Keep track of fluid nodes vs solid nodes
    %Initial Condition Poiseuille flow
    t = 0;
    R = 1;
    tol = 1e-5;
    old_Mag = sqrt(Ux.^2+Uy.^2);
    f = zeros([9,size(Ux)]);
    %Inlet Profile
    Uyin = (6*U0/(H)^2).*x'.*((H)-x');
    Uy(:,1) = Uyin;
    Uxin = zeros(size(X(:,1)));
    for i = 1:length(cx) %Initialize Populations
        CdotU = cx(i).*Ux + cy(i).*Uy;
        UdotU = Ux.^2+Uy.^2;
        f(i,:,:) = rho.*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
    end
    feq = f;
    fnew = f;
    figure(1);
    axis equal
    while t < maxT && R > tol
        t = t+1;
        rho = squeeze(sum(fnew, 1));
        Ux = squeeze(sum(fnew .* reshape(cx, 9, 1, 1), 1)) ./ rho;
        Uy = squeeze(sum(fnew .* reshape(cy, 9, 1, 1), 1)) ./ rho;
        %Collision
        for i = 1:length(cx)
            CdotU = cx(i).*Ux + cy(i).*Uy;
            UdotU = Ux.^2+Uy.^2;
            feq(i,:,:) = rho.*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
        end
        fmirror = 2.*feq - f;
        fnew = (1-beta)*f + beta*fmirror;
        %Streaming
        for i = 1:length(cx)
            fnew(i,:,:) = circshift(fnew(i,:,:),[0,cx(i),cy(i)]);
        end
        %Boundary Conditions Inlet
        rhonew = squeeze(sum(fnew, 1));
        for i = 1:length(cx)
            CdotU = cx(i).*Uxin + cy(i).*Uyin;
            UdotU = Uxin.^2+Uyin.^2;
            fnew(i,:,1) = rho(:,1).*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
        end

        %Boundary Conditions Outlet
        for i= 1:length(cx)
            fnew(i,:,end) = fnew(i,:,end-1);
        end
        %BounceBack on Cyllinder
        for i=1:length(cx)
            fnew(i,boundary) = fnew(opp(i),boundary);
        end
        %Wall BounceBack
        for i=1:length(cx)
            fnew(i,1,:) = f(opp(i),1,:);
            fnew(i,end,:) = f(opp(i),end,:);
        end
        %Macroscopic Values
        f = fnew;
        if t > 2500
            R = sqrt(sum((Uy(:,end) - UyOld(:,end)).^2,'all'))/sqrt(sum(Uy(:,end),'all'));
            % if t > 100 && abs(R-Rold) < 1e-9
            %     break
            % end
        end
        UxOld = Ux;
        UyOld = Uy;
        rhoOld = rho;
        if mod(t,500) == 0
            mag = reshape(sqrt(Ux.^2+Uy.^2),length(x),length(y));
            mag(boundary) = NaN;
            figure(1)
            clf
            imagesc(mag)
            hold on
            quiver(Y,X,Uy,Ux,2,'Color',[0,0,0],'LineWidth',1);
            colorbar
            set(gcf, 'Units', 'normalized', 'Position', [0.0 0.25 1 1]);
            axis equal off;
            title(['Velocity Magnitude with D=',num2str(D)],['t=',num2str(t)])
            p = nsidedpoly(1000, 'Center', [2*D+1,2*D+1], 'Radius', D/2);
            plot(p,'FaceAlpha',1,'FaceColor',[0,0,0]);
            xlim([0,D*5]);
            % figure(2)
            % set(gcf, 'Units', 'normalized', 'Position', [0.75 0.45 0.25 0.45]);
            % title("Recirculation Plot");
            % plot(y,Uy(floor(length(x)/2),:));
            % pause(.1);
            % figure(3)
            % clf;
            % hold on
            % set(gcf, 'Units', 'normalized', 'Position', [0.25 1,0.25 0.45]);
            % plot(Uyin,x,'ro-');
            % plot(Uy(:,end), x, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
            % figure(4);
            % clf;
            % colormap(redblue);
            % axis equal off;
            % set(gcf, 'Units', 'normalized', 'Position', [0.5 0.75,1 1]);
            % imagesc(mag(center(1)-D:center(1)+D,center(1)-D:center(1)+4*D));
            % hold on
            % 
            % colorbar;
            % title(['Velocity Magnitude with D=',num2str(D)],['t=',num2str(t)])
            % tmp1 = center(1)-D:center(1)+D;
            % tmp2 = center(1)-D:center(1)+4*D;
            % p = nsidedpoly(1000, 'Center', [center(1)-tmp1(1)+2,center(1)-tmp2(1)+2], 'Radius', D/2);
            % plot(p,'FaceAlpha',1,'FaceColor',[0,0,0]);
            % axis equal off
        end
    end
    % streamline(Uy,Ux,Y,X,'color',[0,0,0],'LineWidth',0.2)
    Uy_mid = Uy(floor(length(x)/2),:);
    [~,idx] = min(Uy_mid);
    [~,recirculateMax] = min(abs(Uy_mid(idx:end)));
    recirculateEnd = idx+recirculateMax-1;
    recirculationDist(iD) = (recirculateEnd-2.5*D)/D;
end