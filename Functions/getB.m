function [B] = getB(ndof, Mesh, Quad, ip)
% Get B matrix
% Input parameters: ndof = number of degrees of freedom of the element
%                   Quad = quadrature points
%                   ip = number of the current integration point

% shape functions derivatives
dN = zeros(1, ndof);

% Jacobian
J = Mesh.h/2;

% shape functions vector according to the number of DOFs
switch ndof
    case 2
        dN(1,1) = -1/2;
        dN(1,2) = 1/2;

    case 3
        dN(1,1) = Quad.csi(ip,1) - 1/2;
        dN(1,2) = -2*Quad.csi(ip,1);
        dN(1,3) = Quad.csi(ip,1) + 1/2;
end

% B matrix
B = dN./J;

end