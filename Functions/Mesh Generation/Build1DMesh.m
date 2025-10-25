% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function Mesh = Build1DMesh(nsd, ne, L, type, field)
% Build uniform mesh

%% Mesh properties
Mesh.ne = ne; % number of elements
Mesh.L = L; % column length
Mesh.nsd = nsd; % number of spatial dimensions
Mesh.type = type; % type of element
Mesh.field = field; % variable field (u, p, n)

switch type
    case 'L2'
        % number of nodes per element
        Mesh.nne = 2;
        % number of DOFs per element
        Mesh.nDOFe = 2;
    case 'L3'
        % number of nodes per element
        Mesh.nne = 3;
        % number of DOFs per element
        Mesh.nDOFe = 3;
    case 'L4'
        % number of nodes per element
        Mesh.nne = 4;
        % number of DOFs per element
        Mesh.nDOFe = 4;
end

% number of DOFs
Mesh.nDOF = ne * (Mesh.nDOFe - 1) + 1;
% DOFs
Mesh.DOF = (1:Mesh.nDOF).';
% nodal coordinates
Mesh.coords = (0:(L/(Mesh.nDOF-1)):L).';
% total number of nodes
Mesh.nn = length(Mesh.coords);
% connectivity
Mesh.conn = zeros(Mesh.ne, Mesh.nDOFe);
for e = 1:Mesh.ne
    if e == 1
        Mesh.conn(1,:) = 1:Mesh.nDOFe;
    else
        Mesh.conn(e,:) = Mesh.conn(e-1,end):(Mesh.conn(e-1,end) + Mesh.nDOFe -1);
    end
end

% node sets
Mesh.xdofs = 1:Mesh.nDOF;
Mesh.ydofs = zeros(length(Mesh.xdofs),1);

end