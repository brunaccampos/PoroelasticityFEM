clearvars
clear
clc

L = 1;
h = 0.1;
x = (0:h:L);
xd = x./L;

p0 = 1e3;

E = 3e4;
nu = 0.2;
mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));
eta = 1e-3; % fluid viscosity
k = 1e-10; % permeability

t = 0;
dt = h;
tend = 10;

aux1 = 0;
aux2 = 0;

N = 100;

while t <= tend
    td = (lambda + 2*mu)*k * t / (eta * L^2);

    for n = 0:N
        M = pi()*(2*n+1)/2;
        aux1 = aux1 + 2* cos(M*xd)* exp(-M^2 * td)/M^2;
        aux2 = aux2 + 2* sin(M*xd)* exp(-M^2 * td)/M;
    end

    % displacement
    ud = 1 - xd - aux1;
    u = ud*p0*L/(lambda + 2*mu);

    % stress
    sigmad = -1 + aux2;
    sigma = sigmad*p0;

    % pressure
    pd = aux2;
    p = pd*p0;

    t = t + dt;

%     subplot(1,3,1);
%     plot(xd, p); % pressure
%     hold on
%     xlabel('Column depth (m)');
%     ylabel('Pressure (Pa)');
% %     title('Pressure at time t = %d', t);
%     axis([0 1 0 120000]);
%     hold off
% 
%     subplot(1,3,2);
%     plot(xd,sigma); % stress
%     hold on
%     xlabel('Column depth (m)');
%     ylabel('Stress (Pa)');
% %     title('Pressure at time t = %d', t);
%     hold off
% 
%     subplot(1,3,3);
%     plot(xd,u); % displacement
%     hold on
%     xlabel('Column depth (m)');
%     ylabel('Displacement (m)');
% %     title('Pressure at time t = %d', t);
%     hold off
% 
%     pause(0.01)

end

figure
plot(xd,sigma);