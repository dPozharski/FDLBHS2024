% Parameters
icFlag = 1; %Flag for which initial condition to simulate
%0 -> Gaussian hill
%1 -> Hyperbollic Tangent
quadraticTerm = 0;
nx = 500;          % Number of lattice nodes
nt = 5000;        % Number of time steps
nu = 0.05;       % Kinematic viscosity
dx = 1.0;         % Lattice spacing
dt = 1.0;         % Time step
cs = dx/(sqrt(3)*dt); % Lattice speed of sound
c = [-1; 0; 1]; %Velocities
w = [1/6; 2/3; 1/6]; %Weights
beta = (2*(nu/(cs^2*dt) + 1/2))^-1; % LB kinematic viscosity
% Initialize the distribution function
% f = zeros(nx,3); % f(i,j) for i=1:nx and j=1:3 (3 velocity directions)
feq = zeros(nx,3);
rho = zeros(nx,1);
u = zeros(nx,1) + 0.1;
% Initial conditions
rho_0 = 1;
drho = 0.5;
sigma = 1;
x_0 = 250;
switch icFlag
    case 0
        rho = (rho_0 + drho.*exp(-(((1:nx) - x_0)./sigma).^2))';
        maxrhoInit = max(rho);
    case 1
        rho = (rho_0 + .5 .* drho .* (1 + tanh(((1:nx) - x_0)./sigma)))';
        maxrhoInit = max(rho);
end
if quadraticTerm
    % Second order equilibrium (6)
    feq(:,1) = w(1).*rho.*(1+(c(1).*u)./(cs^2) + ((c(1)^2 - cs^2).*u.^2)/(2*cs^4));
    feq(:,2) = w(2).*rho.*(1+(c(2).*u)./(cs^2) + ((c(2)^2 - cs^2).*u.^2)/(2*cs^4));
    feq(:,3) = w(3).*rho.*(1+(c(3).*u)./(cs^2) + ((c(3)^2 - cs^2).*u.^2)/(2*cs^4));
else
    % % First order equilibrium (7)
    feq(:,1) = w(1).*rho.*(1+(c(1).*u(:))./(cs^2));
    feq(:,2) = w(2).*rho.*(1+(c(2).*u(:))./(cs^2));
    feq(:,3) = w(3).*rho.*(1+(c(3).*u(:))./(cs^2));
end
figure(2)
% hold on
f = feq; %Initial condition
fnew = f; %Temp value to hold next timestep

maxrhoInit = max(sum(f,2));
% Time-stepping loop
for t = 1:nt

    rho = sum(f,2); % (eq 4)
    % u = (c(1).*f(:,1) + c(3).*f(:,3))./rho; % (eq 5)

    if quadraticTerm
        % Second order equilibrium (eq 6)
        feq(:,1) = w(1).*rho.*(1+(c(1).*u)./(cs^2) + ((c(1)^2 - cs^2).*u.^2)/(2*cs.^4));
        feq(:,2) = w(2).*rho.*(1+(c(2).*u)./(cs^2) + ((c(2)^2 - cs^2).*u.^2)/(2*cs.^4));
        feq(:,3) = w(3).*rho.*(1+(c(3).*u)./(cs^2) + ((c(3)^2 - cs^2).*u.^2)/(2*cs.^4));
    else
        % % First order equilibrium (eq 7)
        feq(:,1) = w(1).*rho.*(1+(c(1).*u(:))./(cs^2));
        feq(:,2) = w(2).*rho.*(1+(c(2).*u(:))./(cs^2));
        feq(:,3) = w(3).*rho.*(1+(c(3).*u(:))./(cs^2));
    end
    % Collision step
    % fmirror =  2.*feq - f; % (eq 2)
    % f = (1-beta).*f + beta.*fmirror;
    f = f + 2*beta.*(feq-f);
    fnew = f;
    %Streaming step
    fnew(1:end-1,1) = f(2:end,1); % c = -1
    fnew(2:end,3) = f(1:end-1,3); % c = 1

    %Boundary conditions
    switch icFlag
        case 0
            fnew(1,3) = f(end,3);
            fnew(end,1) = f(1,1);
        case 1
            fnew(1,3) = f(1,3);
            fnew(end,1) =  f(end,1);
    end

    f = fnew;

    % Visualization (optional)
    if mod(t, 10) == 0
        mach = abs(u./cs);
        plot(1:nx,rho);
        xlim([0 nx]);
        ylim([.95 maxrhoInit]);
        title('Gaussian Hill Advection Diffusion',['Time Step: ', num2str(t)]);
        pause(0.01);
    end
end