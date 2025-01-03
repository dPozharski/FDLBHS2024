%%%%%%
%Exercise #1 Lattice Boltzman Methods
%Denis Pozharskiy
%25/09/2024
%%%%%%
setGroot();
%Initialization
% Parameters
quadraticTerm = 1; %Flag for the use of the quadratic term

nx = 800;          % Number of lattice nodes
nt = 500;        % Number of time steps
rho0_left = 1.5;  % Density on the left side
rho0_right = 1.0; % Density on the right side
u_left = 0.0;     % Velocity on the left side
u_right = 0.0;    % Velocity on the right side
nu = 1e-5;     % Kinematic viscosity
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
u = zeros(nx,1);
% Initial conditions
rho(:) = rho0_left;
rho(1:nx>(nx/2)) = rho0_right;


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
figure(1)

f = feq; %Initial condition
fnew = f; %Temp value to hold next timestep

% Time-stepping loop
for t = 1:nt

    rho = sum(f,2); % (eq 4)
    u = (c(1).*f(:,1) + c(3).*f(:,3))./rho; % (eq 5)

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
    fnew(1,3) = f(1,1);
    fnew(end,1) = f(end,3);

    f = fnew;

    % Visualization (optional)
    if mod(t, 20) == 0
        mach = u;
        plot(1:nx,rho,1:nx,mach);
        xlim([0 nx]);
        xlabel('r');
        legend('Density','Velocity')
        ylim([0 max(rho)]);
        title('Density and Velocity for 1D Isothermal Shock Tube','Using second order equilibrium');
        pause(0.01);
    end
end
[~,idx] = max(abs(diff(rho(400:end))));
idx = idx + 400;
shockSpeed = ((idx-(nx/2))/(nt-1));
cComp = speedOfSound(rho(nx/2),rho(end),shockSpeed);

function c = speedOfSound(rho2,rho1,vShock)
    gamma = 1.4;
    M = sqrt((2*rho2)/(rho1*((gamma+1)-(rho2/rho1)*(gamma-1)))); %Rankine Hugoniot relations
    c = vShock/M;
end