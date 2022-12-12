function [N] = getN(ndof, Quad, ip)
% Get Shape Functions
% Input parameters: ndof = number of degrees of freedom of the element
%                   Quad = quadrature points
%                   ip = number of the current integration point

% initialize shape function vector
N = zeros(1, ndof);

% shape functions vector according to the number of DOFs
switch ndof
    case 2
        N(1,1) = 0.5*(1 - Quad.csi(ip,1));
        N(1,2) = 0.5*(1 + Quad.csi(ip,1));

    case 3
        N(1,1) = 0.5*(Quad.csi(ip,1)^2 - Quad.csi(ip,1));
        N(1,2) = 1 - Quad.csi(ip,1)^2;
        N(1,3) = 0.5*(Quad.csi(ip,1)^2 + Quad.csi(ip,1));
end

end