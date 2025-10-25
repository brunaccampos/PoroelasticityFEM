% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function Mesh = BuildMesh_structured(nsd, coord0, L, ne, type, field, progress_on)
% Structured mesh generator
% Acknowledgments: Chris Ladubec, Matin Parchei Esfahani

%% Mesh properties
Mesh.type = type;
% variable field (u, p, n)
Mesh.field = field;
% total number of elements
Mesh.ne = prod(ne);
% number of nodes per element (per direction)
nne = NodesPerElement(Mesh.type);
% [nnex;nney] nnex: number of nodes per element in the x-direction
% nney: number of nodes per element in the y-direction
% number of nodes per element
Mesh.nne = prod(nne);
% number of spatial dimensions
Mesh.nsd = nsd;

%% Nodal spatial locations

if ~isfield(Mesh,'x')
    % compute nodal locations
    if progress_on
        disp([num2str(toc),': Computing nodal locations...']);
    end
    dL = L./(ne.*(nne-1));  % node spacing
    x_temp = cell(1,nsd);
    for i = 1:nsd
        % nodal spatial locations
        x_temp{i} = coord0(i):dL(i):(coord0(i)+L(i));
        % total number of nodes per direction
        nnx(i) = length(x_temp{i});
    end

    Mesh.nn = prod(nnx);               % total number of nodes

    Mesh.coords = zeros(Mesh.nn,nsd);
    switch nsd
        case 1
            Mesh.coords = x_temp{1}';   % nodal spatial locations
        case 2
            node_count = 1;
            for nx = 1:nnx(1)
                for ny = 1:nnx(2)
                    Mesh.coords(node_count,1) = x_temp{1}(nx);
                    Mesh.coords(node_count,2) = x_temp{2}(ny);
                    node_count = node_count + 1;
                end
            end
        case 3
            node_count = 1;
            for nx = 1:nnx(1)
                for ny = 1:nnx(2)
                    for nz = 1:nnx(3)
                        Mesh.coords(node_count,1) = x_temp{1}(nx);
                        Mesh.coords(node_count,2) = x_temp{2}(ny);
                        Mesh.coords(node_count,3) = x_temp{3}(nz);
                        node_count = node_count + 1;
                    end
                end
            end
    end
end

%% Nodal matrix

switch Mesh.nsd
    case 1
        % [nodex1 nodex2 ... nodexnn]
        nodes = 1:nnx;
    case 2
        % Matrix of nodes in the x- and y-direction
        % [nodex1y1     nodex2y1    ... nodex_nnxy1;
        % [nodex1y2     nodex2y2    ... nodex_nnxy2;
        % [...;
        % [nodex1y_nny  nodex2y_nny ... nodex_nnx y_nny];
        nodes = zeros(nnx(2),nnx(1));

        node_count = 1;
        for nx = 1:nnx(1)
            for ny = 1:nnx(2)
                nodes(ny,nx) = node_count;
                node_count = node_count + 1;
            end
        end

    case 3
        % 3D Array of nodes in the x-, y- and z-directions
        % [nodex1y1z1     nodex2y1z1    ... nodex_nnxy1z1;
        % [nodex1y2z1     nodex2y2z1    ... nodex_nnxy2z1;
        % [...;
        % [nodex1y_nnyz1  nodex2y_nnyz1 ... nodex_nnx y_nnyz1];
        % ... extended in the third dimension by z2, z3...
        nodes = zeros(nnx(2),nnx(1),nnx(3));
        node_count = 1;
        for nx = 1:nnx(1)
            for ny = 1:nnx(2)
                for nz = 1:nnx(3)
                    nodes(ny,nx,nz) = node_count;
                    node_count = node_count + 1;
                end
            end
        end
end

%% Element connectivity and neighbours
if progress_on
    disp([num2str(toc),': Defining element connectivity...']);
end
% initialize member connectivity
Mesh.conn = zeros(Mesh.ne, Mesh.nne);

% number of neighbors
if nsd == 1
    nneighbors = 2;
else
    nneighbors = 4;
end

% initialize array of neighboring elements
% [e1 neighbor 1, e1 neighbor 2, e1 neighbor 3, e1 neighbor 4 ...;]
% [e2 neighbor 1, e2 neighbor 2, e2 neighbor 3, e2 neighbor 4 ...;]
% [...]
% [e_ne neighbor 1, e_ne neighbor 2, e_ne neighbor 3, e_ne neighbor 4 ...;]

Mesh.eneighbours = zeros(Mesh.ne,nneighbors);

if nsd == 1
    Mesh.L = max(Mesh.coords,[],'all');
    for e = 1:Mesh.ne   % loop through all elements
        % if it is the first element
        if e == 1
            Mesh.conn(1,:) = 1:Mesh.nne;
        else
            Mesh.conn(e,:) = Mesh.conn(e-1,end): ...
                (Mesh.conn(e-1,end) + nne - 1);
        end

        Mesh.eneighbours(e,1) = e - 1;
        Mesh.eneighbours(e,2) = e + 1;
    end
elseif strcmp(Mesh.type,'Q4')
    e = 1; % element counter
    % run through the elements in the x direction
    for nx = 1:ne(1)
        % run through the elements in the y direction
        for ny = 1:ne(2)
            % starting node number for the element
            L1 = [ny,nx];
            L2 = L1 + [0,1];
            L3 = L1 + [1,1];
            L4 = L1 + [1,0];

            Mesh.conn(e,1) = nodes(L1(1),L1(2));
            Mesh.conn(e,2) = nodes(L2(1),L2(2));
            Mesh.conn(e,3) = nodes(L3(1),L3(2));
            Mesh.conn(e,4) = nodes(L4(1),L4(2));

            Mesh.eneighbours(e,1) = e - ne(2);
            Mesh.eneighbours(e,2) = e + 1;
            Mesh.eneighbours(e,3) = e + ne(2);
            Mesh.eneighbours(e,4) = e - 1;

            badleftneighbours = Mesh.eneighbours(e,1) < ny;
            badrightneighbours = Mesh.eneighbours(e,3)> ...
                Mesh.ne - (ne(2) - ny);
            Mesh.eneighbours(e,1) = Mesh.eneighbours(e,1)*(1-badleftneighbours);
            Mesh.eneighbours(e,3) = Mesh.eneighbours(e,3)*(1-badrightneighbours);
            e = e + 1;

        end
        % Set non-existent neighbours to zero (for edge elements)
        % elements in this column
        colel = (nx-1)*ne(2)+1:nx*ne(2);
        badbotneighbours = find(Mesh.eneighbours(colel,4)<(nx-1)*ne(2)+1);
        badtopneighbours = find(Mesh.eneighbours(colel,2)>nx*ne(2));
        Mesh.eneighbours(colel(badbotneighbours),4) = zeros(size(badbotneighbours));
        Mesh.eneighbours(colel(badtopneighbours),2) = zeros(size(badtopneighbours));
    end
elseif strcmp(Mesh.type,'Q9')
    e = 1; % element counter
    % run through the elements in the y-direction
    for ny = 1:ne(2)
        % run through the elements in the x-direction
        for nx = 1:ne(1)
            % starting node
            L1 = 2*(nx-1)*nnx(2) + 2*(ny-1) + 1;
            
            Mesh.conn(e,1) = L1;
            Mesh.conn(e,2) = L1 + 2*nnx(2);
            Mesh.conn(e,3) = L1 + 2 + 2*nnx(2);
            Mesh.conn(e,4) = L1 + 2;
            Mesh.conn(e,5) = L1 + nnx(2);
            Mesh.conn(e,6) = L1 + 1 + 2*nnx(2);
            Mesh.conn(e,7) = L1 + 2 + nnx(2);
            Mesh.conn(e,8) = L1 + 1;
            Mesh.conn(e,9) = L1 + 1 + nnx(2);
            
            Mesh.eneighbours(e,1) = e - ne(1);
            Mesh.eneighbours(e,2) = e + 1;
            Mesh.eneighbours(e,3) = e + ne(1);
            Mesh.eneighbours(e,4) = e - 1;
            
            e = e + 1;

            % Set non-existent neighbours to zero (for edge elements)
            Mesh.eneighbours(Mesh.eneighbours <= 0) = ...
                zeros(size(find(Mesh.eneighbours <= 0)));

            Mesh.eneighbours(Mesh.eneighbours > Mesh.ne) = ...
                zeros(size(find(Mesh.eneighbours > Mesh.ne)));

        end
        % Set non-existent neighbours to zero (for edge elements)
        Mesh.eneighbours(Mesh.eneighbours <= 0) = ...
            zeros(size(find(Mesh.eneighbours <= 0)));

        Mesh.eneighbours(Mesh.eneighbours > Mesh.ne) = ...
            zeros(size(find(Mesh.eneighbours > Mesh.ne)));

    end
elseif strcmp(Mesh.type,'B8')
    e = 1;
    for nx = 1:ne(1)
        for ny = 1:ne(2)
            for nz = 1:ne(3)
                L1 = [ny,nx,nz];
                L2 = L1 + [0,1,0];
                L3 = L1 + [1,1,0];
                L4 = L1 + [1,0,0];
                L5 = L1 + [0,0,1];
                L6 = L1 + [0,1,1];
                L7 = L1 + [1,1,1];
                L8 = L1 + [1,0,1];

                Mesh.conn(e,1) = nodes(L1(1),L1(2),L1(3));
                Mesh.conn(e,2) = nodes(L2(1),L2(2),L2(3));
                Mesh.conn(e,3) = nodes(L3(1),L3(2),L3(3));
                Mesh.conn(e,4) = nodes(L4(1),L4(2),L4(3));
                Mesh.conn(e,5) = nodes(L5(1),L5(2),L5(3));
                Mesh.conn(e,6) = nodes(L6(1),L6(2),L6(3));
                Mesh.conn(e,7) = nodes(L7(1),L7(2),L7(3));
                Mesh.conn(e,8) = nodes(L8(1),L8(2),L8(3));

                e = e + 1;
            end
        end
    end
end

%% Nodal DOFs
Mesh = NodeDOFs(Mesh);

%% Nodal Connectivity
% list of elements connected to each node
Mesh.nodeconn = NodalConn(Mesh);

%% Node sets
Mesh = NodeSets(Mesh);

if progress_on
    disp([num2str(toc),': Done generating mesh...']);
end
end