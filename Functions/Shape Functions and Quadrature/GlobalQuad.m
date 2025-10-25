% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Quad] = GlobalQuad(Mesh, Control)
% Return quadrature points (coordinates and weights), shape functions and 
% its derivatives
% Reference: Fish, J., & Belytschko, T. (2007). A first course in finite
% elements. Chichester; Hoboken, NJ: John Wiley & Sons.
% Acknowledgments: Matin Parchei Esfahani

% number of spatial dimensions
nsd = Mesh.nsd;

% quadrature order
switch Mesh.field
    case 'u'
        nq = Control.nqU;
    otherwise
        nq = Control.nqP;
end

% quadrature points
Quad = getIPs(nsd, nq, Mesh.type);

% number of quadrature points
Quad.nq = length(Quad.w);

% cell array of shape functions at each quadrature point
Quad.Nq = cell(Quad.nq,1);

% cell array of derivative of shape functions at each quadrature point 
Quad.dNdxiq = cell(Quad.nq,1);
Quad.Nv = cell(Quad.nq,1);

for q = 1:Quad.nq
    % quadrature point in parent coordinate
    coord = Quad.p(q,:);       
    % N: shape function evaluated at quadrature point
    % dNdxi: derivative of shape function wrt parent coordinate, 
    %          xi, evaluated at quadrature point
    [Quad.Nq{q}, Quad.dNdxiq{q}] = lagrange_basis(Mesh, coord);  
    Quad.Nv{q} = getNVoigt(Mesh, Quad.Nq{q});
end
