% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [dN] = getdN(Mesh, Quad, ip)
% Get Derivatives of Shape Functions

% number of nodes per element
nne = Mesh.nne;
% element type
type = Mesh.type;

% initialize shape function vector
dN = zeros(Mesh.nsd, nne);

% shape functions vector according to the number of DOFs
switch type
    case 'L2'
        % 1D linear element
        dN(1,1) = -1/2;
        dN(1,2) = 1/2;
    case 'L3'
        % 1D quadratic element
        dN(1,1) = 0.5*(2*Quad.p(ip) - 1);
        dN(1,2) = -2*Quad.p(ip);
        dN(1,3) = 0.5*(2*Quad.p(ip) + 1);
    case 'L4'
        % 1D cubic element
        dN(1,1) = (9*Quad.p(ip))/8 - (27*Quad.p(ip)^2)/16 + 1/16;
        dN(1,2) = (81*Quad.p(ip)^2)/16 - (9*Quad.p(ip))/8 - 27/16;
        dN(1,3) = 27/16 - (81*Quad.p(ip)^2)/16 - (9*Quad.p(ip))/8;
        dN(1,4) = (9*Quad.p(ip))/8 + (27*Quad.p(ip)^2)/16 - 1/16;
    case 'T3'
        % 2D linear element
        % csi derivatives
        dN(1,1) = -1;
        dN(1,2) = 1;
        dN(1,3) = 0;
        % eta derivatives
        dN(2,1) = -1;
        dN(2,2) = 0;
        dN(2,3) = 1;
    case 'T6'
        % 2D quadratic element
        csi = Quad.p(:,1);
        eta = Quad.p(:,2);
        % csi derivatives
        dN(1,1) = -3 + 4 * eta(ip) + 4 * csi(ip);
        dN(1,2) = 4 * csi(ip) - 1;
        dN(1,3) = 0;
        dN(1,4) = 4 - 8 * csi(ip) - 4 * eta(ip);
        dN(1,5) = 4 * eta(ip);
        dN(1,6) = -4 * eta(ip);
        % eta derivatives
        dN(2,1) = -3 + 4 * csi(ip) + 4 * eta(ip);
        dN(2,2) = 0;
        dN(2,3) = 4 * eta(ip) - 1;
        dN(2,4) = -4 * csi(ip);
        dN(2,5) = 4 * csi(ip);
        dN(2,6) = 4 - 4 * csi(ip) - 8 * eta(ip);
    case 'Q4'
        % 2D linear element
        csi = Quad.p(:,1);
        eta = Quad.p(:,2);
        % csi derivatives
        dN(1,1) = -0.25*(1-eta(ip));
        dN(1,2) = 0.25*(1-eta(ip));
        dN(1,3) = 0.25*(1+eta(ip));
        dN(1,4) = -0.25*(1+eta(ip));
        % eta derivatives
        dN(2,1) = -0.25*(1-csi(ip));
        dN(2,2) = -0.25*(1+csi(ip));
        dN(2,3) = 0.25*(1+csi(ip));
        dN(2,4) = 0.25*(1-csi(ip));
    case 'Q9'
        % 2D quadratic element
        csi = Quad.p(:,1);
        eta = Quad.p(:,2);
        % csi derivatives
        dN(1,1) = 0.25*(eta(ip)^2 - eta(ip))*(2*csi(ip) - 1);
        dN(1,2) = 0.25*(eta(ip)^2 - eta(ip))*(2*csi(ip) + 1);
        dN(1,3) = 0.25*(eta(ip)^2 + eta(ip))*(2*csi(ip) + 1);
        dN(1,4) = 0.25*(eta(ip)^2 + eta(ip))*(2*csi(ip) - 1);
        dN(1,5) = 0.5*(eta(ip)^2 - eta(ip))*(-2*csi(ip));
        dN(1,6) = 0.5*(1 - eta(ip)^2)*(2*csi(ip) + 1);
        dN(1,7) = 0.5*(eta(ip)^2 + eta(ip))*(-2*csi(ip));
        dN(1,8) = 0.5*(1 - eta(ip)^2)*(2*csi(ip) - 1);
        dN(1,9) = (1 - eta(ip)^2)*(-2*csi(ip));
        % eta derivatives
        dN(2,1) = 0.25*(csi(ip)^2 - csi(ip))*(2*eta(ip) - 1);
        dN(2,2) = 0.25*(csi(ip)^2 + csi(ip))*(2*eta(ip) - 1);
        dN(2,3) = 0.25*(csi(ip)^2 + csi(ip))*(2*eta(ip) + 1);
        dN(2,4) = 0.25*(csi(ip)^2 - csi(ip))*(2*eta(ip) + 1);
        dN(2,5) = 0.5*(1 - csi(ip)^2)*(2*eta(ip) - 1);
        dN(2,6) = 0.5*(csi(ip)^2 + csi(ip))*(-2*eta(ip));
        dN(2,7) = 0.5*(1 - csi(ip)^2)*(2*eta(ip) + 1);
        dN(2,8) = 0.5*(csi(ip)^2 - csi(ip))*(-2*eta(ip));
        dN(2,9) = (1-csi(ip)^2)*(-2*eta(ip));
end
end