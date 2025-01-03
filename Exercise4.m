%%%%%%
%Exercise #4 Lattice Boltzman Methods
%Denis Pozharskiy
%23/12/2024
%%%%%%
clear;
close all;
setGroot();
%Initialization and Constants
alpha = 80;
delta = 0.05;
L = 512;
[x,y] = ndgrid(0:L-1,0:L-1);
dr = 1;
dt = 1;
Q = 9;
T0 = 1/3;
rho0 = 1;
P0 = rho0/3;
cs = sqrt(T0);
U0 = linspace(0.05,0.5,50);
nu = linspace(0.0001,0.05,50);
[U0,Nu] = ndgrid(U0,nu);
stableMapBGK = zeros(size(U0));
stableMapKBC = zeros(size(U0));
load('KBCSym.mat')
KBC.feqH = matlabFunction(feqH);
KBC.feqK = matlabFunction(feqK);
KBC.feqS = matlabFunction(feqS);
KBC.fH = matlabFunction(fH);
KBC.fK = matlabFunction(fK);
KBC.fS = matlabFunction(fS);
KBC.gamma = matlabFunction(gamma);
numCases = length(U0(:));
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,numCases) '\n\n']);
parfor iC = 1:length(U0(:))
    u0 = U0(iC);
    nu = Nu(iC);
    u0 = 0.1;
    nu = 0.0008;
    ux = u0 * tanh(alpha*(.25 - abs(y/L - 0.5)));
    uy = u0 * delta * sin(2*pi * (x/L + 0.25));
    P = P0 + 0.*x;
    tprime = 5;
    t = tprime*L/u0;
    beta = (2*nu/T0 + 1)^-1;
    rho = P /cs^2;
    KEi = sum(rho.*(ux.^2+uy.^2),'all');
    tc = L/u0;
    [rho,ux,uy] = LBMBGK(ux,uy,rho,t,beta,x,y,tc,nu,u0);
    KE = sum(rho.*(ux.^2+uy.^2),'all');
    if KE <= KEi
        stableMapBGK(iC) = -1;
    end
    ux = u0 * tanh(alpha*(.25 - abs(y/L - 0.5)));
    uy = u0 * delta * sin(2*pi * (x/L + 0.25));
    P = P0 + 0.*x;
    rho = P /cs^2;
    [rho,ux,uy] = LBMKBC(ux,uy,rho,t,beta,x,y,KBC);
    KE = sum(rho.*(ux.^2+uy.^2),'all');
    if KE <= KEi
        stableMapKBC(iC) = -1;
    end
    fprintf('\b|\n');
end




function [rho,Ux,Uy] = LBMBGK(Ux0,Uy0,rho,tMax,beta,x,y,tc,nu,U0)
cx = [0 0 1 1 1 0 -1 -1 -1];
cy = [0 1 1 0 -1 -1 -1 0 1];
w = [4/9 1/9 1/36 1/9 1/36 1/9 1/36 1/9 1/36];
T0 = 1/3;
cs = sqrt(T0);
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
    if t == -1
        figure(1)
        colormap(redblue);
        axis off
        clf
        uMag = sqrt(Ux.^2 + Uy.^2);
        vorticity = curl(Ux,Uy);
        imagesc(flip(vorticity'))
        title(['Vorticity at t=t_c with L = ',num2str(length(x(1,:)))],['nu = ',num2str(nu),' U_0 = ',num2str(U0)])
        %hold on
        %quiver(x,y,Ux,Uy);
        pause(.001);
    end
end
end


function [rho,Ux,Uy] = LBMKBC(Ux0,Uy0,rho,tMax,beta,x,y,KBC)
cx = [0 0 1 1 1 0 -1 -1 -1];
cy = [0 1 1 0 -1 -1 -1 0 1];
w = [4/9 1/9 1/36 1/9 1/36 1/9 1/36 1/9 1/36];
T0 = 1/3;
L = length(x);
cs = sqrt(T0);
f = zeros([size(Ux0),9]);
feq = f;
Ux = Ux0;
Uy = Uy0;
feqK = permute(reshape(cell2mat(arrayfun(KBC.feqK,rho,Ux,Uy,'UniformOutput',false)),L,9,L),[1,3,2]);
feqH = permute(reshape(cell2mat(arrayfun(KBC.feqH,rho,Ux,Uy,'UniformOutput',false)),L,9,L),[1,3,2]);
feqS = permute(reshape(cell2mat(arrayfun(KBC.feqS,rho,Ux,Uy,'UniformOutput',false)),L,9,L),[1,3,2]);
f = feqK+feqH+feqS;
for t = 1:tMax
    %Collision
    if t ~= 1
    
    f1 = f(:,:,1);
    f2 = f(:,:,2);
    f3 = f(:,:,3);
    f4 = f(:,:,4);
    f5 = f(:,:,5);
    f6 = f(:,:,6);
    f7 = f(:,:,7);
    f8 = f(:,:,8);
    f9 = f(:,:,9);
    %feqK = permute(reshape(cell2mat(arrayfun(KBC.feqK,rho,Ux,Uy,'UniformOutput',false)),L,9,L),[1,3,2]);
    %feqH = permute(reshape(cell2mat(arrayfun(KBC.feqH,rho,Ux,Uy,'UniformOutput',false)),L,9,L),[1,3,2]);
    feqH = reshape(KBC.feqH(rho,Ux,Uy),L,L,9);
    %feqS = permute(reshape(cell2mat(arrayfun(KBC.feqS,rho,Ux,Uy,'UniformOutput',false)),L,9,L),[1,3,2]);
    feqS = reshape(KBC.feqS(rho,Ux,Uy),L,L,9);;
    fK = permute(reshape(cell2mat(arrayfun(KBC.fK,f1,f2,f3,f4,f5,f6,f7,f8,f9,'UniformOutput',false)),L,9,L),[1,3,2]);
    %fH = permute(reshape(cell2mat(arrayfun(KBC.fH,f3,f5,f7,f9,'UniformOutput',false)),L,9,L),[1,3,2]);
    fH = reshape(KBC.fH(f3,f5,f7,f9),L,L,9);
    %fS = permute(reshape(cell2mat(arrayfun(KBC.fS,f2,f3,f4,f5,f6,f7,f8,f9,'UniformOutput',false)),L,9,L),[1,3,2]);
    fS = reshape(KBC.fS(f2,f3,f4,f5,f6,f7,f8,f9),L,L,9);
    % gamma = permute(reshape(cell2mat(arrayfun(KBC.gamma,beta.*ones(size(f1)),f2,f3,f4,f5,f6,f7,f8,f9,rho,Ux,Uy,'UniformOutput',false)),L,1,L),[1,3,2]);
    gamma = KBC.gamma(beta.*ones(size(f1)),f2,f3,f4,f5,f6,f7,f8,f9,rho,Ux,Uy);
    fmirror = fK + (2.*feqS - fS) + ((1-gamma).*fH + gamma.*feqH);
    f = (1-beta)*f + beta*fmirror;
    end
    %Streaming
    for i = 1:length(cx)
        f(:,:,i) = circshift(f(:,:,i),[cx(i),cy(i)]);
    end
    rho = sum(f, 3);
    Ux = sum(f .* reshape(cx, 1, 1, 9), 3) ./ rho;
    Uy = sum(f .* reshape(cy, 1, 1, 9), 3) ./ rho;
    if mod(t,100) == -1
        figure(1)
        clf
        uMag = sqrt(Ux.^2 + Uy.^2);
        vorticity = curl(Ux,Uy);
        surface(x,y,vorticity,'EdgeColor','none');
        %hold on
        %quiver(x,y,Ux,Uy);
        pause(.001);
    end
end
end