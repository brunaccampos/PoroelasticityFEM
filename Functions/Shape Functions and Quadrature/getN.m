function [N] = getN(Mesh, Quad, ip)
% Get Shape Functions
% Input parameters: Mesh
%                   Quad = quadrature points
%                   ip = number of the current integration point

% number of nodes per element
nne = Mesh.nne;
% element type
type = Mesh.type;

% initialize shape function vector
N = zeros(1, nne);

% shape functions vector according to the number of DOFs
switch type
    case 'L2'
        % 1D linear element
        N(1,1) = 0.5*(1 - Quad.p(ip,1));
        N(1,2) = 0.5*(1 + Quad.p(ip,1));
    case 'L3'
        % 1D quadratic element
        N(1,1) = 0.5*(Quad.p(ip,1)^2 - Quad.p(ip,1));
        N(1,2) = 1 - Quad.p(ip,1)^2;
        N(1,3) = 0.5*(Quad.p(ip,1)^2 + Quad.p(ip,1));
    case 'T3'
        % 2D linear element
        csi = Quad.p(:,1);
        eta = Quad.p(:,2);
        N(1,1) = 1 - csi(ip) - eta(ip);
        N(1,2) = csi(ip);
        N(1,3) = eta(ip);
    case 'T6'
        % 2D quadratic element
        csi = Quad.p(:,1);
        eta = Quad.p(:,2);
        N(1,1) = 1 - 3*(csi(ip) + eta(ip)) + 4*csi(ip)*eta(ip) + 2*(csi(ip)^2+eta(ip)^2);
        N(1,2) = csi(ip) * (2*csi(ip) - 1);
        N(1,3) = eta(ip) * (2*eta(ip) - 1);
        N(1,4) = 4 * csi(ip) * (1 - csi(ip) - eta(ip));
        N(1,5) = 4 * csi(ip) * eta(ip);
        N(1,6) = 4 * eta(ip) * (1 - csi(ip) - eta(ip));
    case 'Q4'
        % 2D linear element
        csi = Quad.p(:,1);
        eta = Quad.p(:,2);
        N(1,1) = 0.25*(1 - csi(ip))*(1 - eta(ip));
        N(1,2) = 0.25*(1 + csi(ip))*(1 - eta(ip));
        N(1,3) = 0.25*(1 + csi(ip))*(1 + eta(ip));
        N(1,4) = 0.25*(1 - csi(ip))*(1 + eta(ip));
    case 'Q9'
        % 2D quadratic element
        csi = Quad.p(:,1);
        eta = Quad.p(:,2);
        N(1,1) = 0.25*(csi(ip)^2 - csi(ip))*(eta(ip)^2 - eta(ip));
        N(1,2) = 0.25*(csi(ip)^2 + csi(ip))*(eta(ip)^2 - eta(ip));
        N(1,3) = 0.25*(csi(ip)^2 + csi(ip))*(eta(ip)^2 + eta(ip));
        N(1,4) = 0.25*(csi(ip)^2 - csi(ip))*(eta(ip)^2 + eta(ip));
        N(1,5) = 0.5*(eta(ip)^2 - eta(ip))*(1 - csi(ip)^2);
        N(1,6) = 0.5*(csi(ip)^2 + csi(ip))*(1 - eta(ip)^2);
        N(1,7) = 0.5*(eta(ip)^2 + eta(ip))*(1 - csi(ip)^2);
        N(1,8) = 0.5*(csi(ip)^2 - csi(ip))*(1 - eta(ip)^2);
        N(1,9) = (1-csi(ip)^2)*(1-eta(ip)^2);
end
end