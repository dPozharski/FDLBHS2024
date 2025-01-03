%%%%%%
%Exercise #5 Lattice Boltzman Methods
%Denis Pozharskiy
%23/12/2024
%%%%%%
clear;
close all;
setGroot();
%Initialization and Constants
dt =  1;
dr = 1;
Re = 50;
Hs =[32,64,96];
L = 10;
U0 = 0.1;
Q = 9;
T0 = 1/3;
rho0 = 1;
P0 = rho0/3;
cs = sqrt(T0);
Errors = [0,0,0];
for iH = 1:length(Hs)
    H = Hs(iH);
    nu = U0*H/Re;
    Fx = 12*rho0*nu*U0/H^2;
    [x,y] = ndgrid(0:H-1,0:L-1);
    Ux = zeros(size(x));
    Uy = zeros(size(y));
    P = P0 + 0.*x;
    tprime = 150;
    t = tprime*L/U0;
    beta = (2*nu/T0 + 1)^-1;
    rho = P /cs^2;
    xs = linspace(0,H-1,H);
    UxA = (Fx/(2*rho0*nu)).*xs.*((H-1)-xs);
    [rho,Ux,Uy,Error] = LBMBGKForcing(Ux,Uy,rho, 5000000,beta,x,y,0,Fx,UxA,xs);
    Errors(iH) = Error;
end

E = log10(Errors);
drL = log10(Hs);
figure(1)
clf;
plot(drL,E);
title("Error From Analytical Poiseuille Flow");
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
xlabel('$\log(\frac{dH}{H})$','interpreter','latex','FontWeight','bold')
ylabel('$\log(E)$','interpreter','latex','FontWeight','bold')
grid on
function [rho,Ux,Uy,Error] = LBMBGKForcing(Ux0,Uy0,rho,tMax,beta,x,y,Fx,Fy,VelProfile,xs)
cx = [0 0 1 1 1 0 -1 -1 -1];
cy = [0 1 1 0 -1 -1 -1 0 1];
w = [4/9 1/9 1/36 1/9 1/36 1/9 1/36 1/9 1/36];

T0 = 1/3;
cs = sqrt(T0);
f = zeros([size(Ux0),9]);
feq = f;
feqForced = f;
Ux = Ux0;
Uy = Uy0;
for i = 1:length(cx)
    CdotU = cx(i).*Ux + cy(i).*Uy;
    UdotU = Ux.^2+Uy.^2;
    f(:,:,i) = rho.*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
end
t = 0;
r = 1;
while t <tMax && r > 1e-9
    t = t+1;
    %Collision
    for i = 1:length(cx)
        CdotU = cx(i).*Ux + cy(i).*Uy;
        UdotU = Ux.^2+Uy.^2;
        feq(:,:,i) = rho.*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
    end
    fmirror = 2.*feq - f;
    forcedUx = Ux + Fx;
    forcedUy = Uy + Fy;
    for i = 1:length(cx)
        CdotU = cx(i).*forcedUx + cy(i).*forcedUy;
        UdotU = forcedUx.^2+forcedUy.^2;
        feqForced(:,:,i) = rho.*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
    end
    F = feqForced-feq;
    f = (1-beta)*f + beta*fmirror + F;
    %Streaming
    for i = 1:length(cx)
        f(:,:,i) = circshift(f(:,:,i),[cx(i),cy(i)]);
    end
    %Boundary Conditions
    f(1,:,[5,4,3]) = f(1,:,[9,8,7]);
    f(end,:,[9,8,7]) = f(end,:,[5,4,3]);
    rho = sum(f, 3);
    Ux = sum(f .* reshape(cx, 1, 1, 9), 3) ./ rho;
    Uy = sum(f .* reshape(cy, 1, 1, 9), 3) ./ rho;
    ux_real = Uy + Fx./(2*rho);
    if mod(t,50) == 0
        figure(1)
        clf
        set(gcf, 'Units', 'normalized', 'Position', [0.0 0.45 0.25 0.45]);
        hold on
        plot(VelProfile,xs,'ro-','Marker','+','LineStyle','--','LineWidth',4);
        plot(ux_real(:,length(y(1,:)) / 2), x(:,1), 'bo-', 'LineWidth', 2, 'MarkerSize', 6);

    end
    if t > 1
        r = sqrt(sum((ux_4real - ux_real_old).^2,'all'))/sqrt(sum(ux_real,'all'));
    end
    ux_real_old = ux_real;
end
Error = sqrt(sum((ux_real(:,length(y(1,:)) / 2)' - VelProfile).^2))/sqrt(sum(VelProfile.^2));
end