clear all;
close all;
setGroot();
Ls = [64,96,128];
errors = zeros(size(Ls));
ErrorOverTime = cell(size(Ls));
for iL = 1:length(Ls)
    %IC
    L = Ls(iL);
    dt = 1;
    dr = 1;
    rho_0 = 1;
    P_0 = rho_0/3;
    U_0 = 0.01;
    Re = 100;
    nu = U_0*L/Re;
    T0 = 1/3;
    cs = sqrt(T0);
    tprime = 2.5;
    t = tprime*L/U_0;
    beta = (2*nu/T0 + 1)^-1;

    %Grid
    [x,y] = ndgrid(0:L-1,0:L-1);
    ux = -U_0.*cos(2.*pi.*x./L).*sin(2.*pi.*y./L);
    uy = U_0.*sin(2.*pi.*x./L).*cos(2.*pi.*y./L);
    P = P_0 - ((rho_0*U_0^2)/4).*(cos(4.*pi.*x./L)+cos(4.*pi.*y./L));
    rho = P/cs^2;

    [Ux,Uy,E] = LBMBGK(ux,uy,rho,t,beta,x,y,nu,U_0);

    UxA = (ux./U_0).*exp((-8*tprime*pi^2)/Re);
    UyA = (uy./U_0).*exp((-8*tprime*pi^2)/Re);

    UxNd = Ux./U_0;
    UyNd = Uy./U_0;
    UxANd = UxA./U_0;
    UyANd = UyA./U_0;
    ErrorOverTime{iL} = E;

    Error = L2Error(UxNd,UyNd,UxA,UyA);
    errors(iL) = Error;
end

xs = log10(dr./Ls);
ys = log10(errors);

[P,S] = polyfit(xs,ys,1);

% loglog(xs,ys);
figure(1)
set(gcf, 'Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5]);
grid on
hold on
plot(ErrorOverTime{1})
plot(ErrorOverTime{2})
plot(ErrorOverTime{3})
legend("L=64","L=96","L=128",Location='northeast');
xlabel("Time Step");
ylabel("L2 Error");
title(['Convergence With Acoustic Scaling (U_0 = ',num2str(U_0),')']);
figure(3)
E = log10([ErrorOverTime{1}(end),ErrorOverTime{2}(end),ErrorOverTime{3}(end)]);
drL = log10(Ls);
plot(drL,E);
%% Part 2
Ls = [64,96,128];
errors = zeros(size(Ls));
ErrorOverTime = cell(size(Ls));
for iL = 1:length(Ls)
    %IC
    L = Ls(iL);
    dt = 1;
    dr = 1;
    rho_0 = 1;
    P_0 = rho_0/3;
    nu = 0.064;
    Re = 100;
    U_0 = Re*nu/L;
    T0 = 1/3;
    cs = sqrt(T0);
    tprime = 2.5;
    t = tprime*L/U_0;
    beta = (2*nu/T0 + 1)^-1;

    %Grid
    [x,y] = ndgrid(0:L-1,0:L-1);
    ux = -U_0.*cos(2.*pi.*x./L).*sin(2.*pi.*y./L);
    uy = U_0.*sin(2.*pi.*x./L).*cos(2.*pi.*y./L);
    P = P_0 - ((rho_0*U_0^2)/4).*(cos(4.*pi.*x./L)+cos(4.*pi.*y./L));
    rho = P/cs^2;

    [Ux,Uy,E] = LBMBGK(ux,uy,rho,t,beta,x,y,nu,U_0);

    UxA = (ux./U_0).*exp((-8*tprime*pi^2)/Re);
    UyA = (uy./U_0).*exp((-8*tprime*pi^2)/Re);

    UxNd = Ux./U_0;
    UyNd = Uy./U_0;
    UxANd = UxA./U_0;
    UyANd = UyA./U_0;
    ErrorOverTime{iL} = E;

    Error = L2Error(UxNd,UyNd,UxA,UyA);
    errors(iL) = Error;
end
figure(2)
set(gcf, 'Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5]);
hold on
grid on
plot(ErrorOverTime{1})
plot(ErrorOverTime{2})
plot(ErrorOverTime{3})
legend("L=64","L=96","L=128",Location='northeast');
xlabel("Time Step");
ylabel("L2 Error");
title("Convergence With Diffusive Scaling (\nu = 0.064)");
figure(3)
set(gcf, 'Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5]);
hold on;
E = log10([ErrorOverTime{1}(end),ErrorOverTime{2}(end),ErrorOverTime{3}(end)]);
drL = log10(Ls);
plot(drL,E);
title("Order of Convergence");
legend('Acoustic Scaling','Diffusive Scaling');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
grid on
%2D Lattice Boltzman using BGK collision operator
function [Ux,Uy,E] = LBMBGK(Ux0,Uy0,rho,tMax,beta,x,y,nu,U0)
cx = [-1,0,1,-1,0,1,-1,0,1];
cy = [-1,-1,-1,0,0,0,1,1,1];
w = [1/36,1/9,1/36,1/9,4/9,1/9,1/36,1/9,1/36];
T0 = 1/3;
cs = sqrt(T0);
L = length(Ux0);
E = zeros(tMax,1);
f = zeros([size(Ux0),9]);
feq = f;
Ux = Ux0;
Uy = Uy0;
for i = 1:length(cx)
    CdotU = cx(i).*Ux + cy(i).*Uy;
    UdotU = Ux.^2+Uy.^2;
    f(:,:,i) = rho.*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
end
for t = 1:tMax
    %Collision
    for i = 1:length(cx)
        CdotU = cx(i).*Ux + cy(i).*Uy;
        UdotU = Ux.^2+Uy.^2;
        feq(:,:,i) = rho.*w(i).*(1 + CdotU/cs^2 + ((CdotU.*CdotU) - cs^2.*(UdotU))/(2*cs^4));
    end
    fmirror = 2.*feq - f;
    f = (1-beta)*f + beta*fmirror;
    %Streaming
    for i = 1:length(cx)
        f(:,:,i) = circshift(f(:,:,i),[cx(i),cy(i)]);
    end
    rho = sum(f, 3);
    Ux = sum(f .* reshape(cx, 1, 1, 9), 3) ./ rho;
    Uy = sum(f .* reshape(cy, 1, 1, 9), 3) ./ rho;
    if true
        UxA = Ux0.*exp((-8*nu*t*pi^2)/L^2);
        UyA = Uy0.*exp((-8*nu*t*pi^2)/L^2);
        E(t) = L2Error(Ux,Uy,UxA,UyA);
    end
    if mod(t,1) == 1
        figure(1)
        clf
        uMag = sqrt(Ux.^2 + Uy.^2);
        uMagA = sqrt(UxA.^2+UyA.^2);
        surface(x,y,abs(uMag-uMagA),'EdgeColor','none');
        %hold on
        %quiver(x,y,Ux,Uy);
        pause(.001);
    end
end
end
function [Error] = L2Error(Ux,Uy,UxA,UyA)
Error = sqrt((sum((Ux-UxA).^2 + (Uy-UyA).^2,'all')))/sqrt(sum(UxA.^2 + UyA.^2,'all'));
end