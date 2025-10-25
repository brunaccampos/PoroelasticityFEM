% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function Mesh = NodeSets(Mesh)
%   Updates the structure array containing mesh information with relevant 
%   node sets

switch Mesh.nsd
    case 1
        Mesh.left_nodes = find(Mesh.coords(:,1)==min(Mesh.coords(:,1)));
        Mesh.right_nodes = find(Mesh.coords(:,1)==max(Mesh.coords(:,1)));

        Mesh.left_dof = Mesh.DOF(Mesh.left_nodes);
        Mesh.right_dof = Mesh.DOF(Mesh.right_nodes);

        Mesh.xdofs = 1:Mesh.nDOF;
        Mesh.ydofs = zeros(length(Mesh.xdofs),1);
        Mesh.zdofs = [];
        
        % corner nodes (structured: nodes numbered sequentially)
        Mesh.conn_corner_nodes = Mesh.conn(:,[1,end]);
        Mesh.corner_nodes = unique(Mesh.conn_corner_nodes);

    case 2
        switch Mesh.field
            case 'u'
                % Nodes and dofs on the left face of the domain
                % (left and right are defined along the x-direction)
                Mesh.left_nodes = find(Mesh.coords(:,1)==min(Mesh.coords(:,1)));
                Mesh.left_dof = Mesh.DOF(Mesh.left_nodes,:);
                Mesh.left_dofx = Mesh.left_dof(:,1);
                Mesh.left_dofy = Mesh.left_dof(:,2);
                Mesh.left_dof = reshape(Mesh.left_dof',numel(Mesh.left_dof),1);

                % Nodes and dofs on the right face of the domain
                Mesh.right_nodes = find(Mesh.coords(:,1)==max(Mesh.coords(:,1)));
                Mesh.right_dof = Mesh.DOF(Mesh.right_nodes,:);
                Mesh.right_dofx = Mesh.right_dof(:,1);
                Mesh.right_dofy = Mesh.right_dof(:,2);
                Mesh.right_dof = reshape(Mesh.right_dof',numel(Mesh.right_dof),1);

                % Nodes and dofs on the top face of the domain
                % (top and bottom are defined along the y-direction)
                Mesh.top_nodes = find(Mesh.coords(:,2)==max(Mesh.coords(:,2)));
                Mesh.top_dof = Mesh.DOF(Mesh.top_nodes,:);
                Mesh.top_dofx = Mesh.top_dof(:,1);
                Mesh.top_dofy = Mesh.top_dof(:,2);
                Mesh.top_dof = reshape(Mesh.top_dof',numel(Mesh.top_dof),1);

                % Nodes and dofs on the bottom face of the domain
                Mesh.bottom_nodes = find(Mesh.coords(:,2)==min(Mesh.coords(:,2)));
                Mesh.bottom_dof =  Mesh.DOF(Mesh.bottom_nodes,:);
                Mesh.bottom_dofx = Mesh.bottom_dof(:,1);
                Mesh.bottom_dofy = Mesh.bottom_dof(:,2);
                Mesh.bottom_dof = reshape(Mesh.bottom_dof',numel(Mesh.bottom_dof),1);

                % DOFs in the x- and y- directions
                Mesh.xdofs = 1:2:Mesh.nDOF;
                Mesh.ydofs = 2:2:Mesh.nDOF;
                Mesh.zdofs = [];

            case {'p', 'n'}
                % Nodes and dofs on the left face of the domain
                % (left and right are defined along the x-direction)
                Mesh.left_nodes = find(Mesh.coords(:,1)==min(Mesh.coords(:,1)));
                Mesh.left_dof = Mesh.DOF(Mesh.left_nodes);

                % Nodes and dofs on the right face of the domain
                Mesh.right_nodes = find(Mesh.coords(:,1)==max(Mesh.coords(:,1)));
                Mesh.right_dof = Mesh.DOF(Mesh.right_nodes);

                % Nodes and dofs on the top face of the domain
                % (top and bottom are defined along the y-direction)
                Mesh.top_nodes = find(Mesh.coords(:,2)==max(Mesh.coords(:,2)));
                Mesh.top_dof = Mesh.DOF(Mesh.top_nodes);

                % Nodes and dofs on the bottom face of the domain
                Mesh.bottom_nodes = find(Mesh.coords(:,2)==min(Mesh.coords(:,2)));
                Mesh.bottom_dof =  Mesh.DOF(Mesh.bottom_nodes);

                % DOFs in the x- and y- directions
                Mesh.xdofs = 1:Mesh.nDOF;
        end
        
        % corner nodes (GMSH/structured: first 4 columns contains corner nodes)
        Mesh.conn_corner_nodes = Mesh.conn(:,1:4);
        Mesh.corner_nodes = unique(Mesh.conn_corner_nodes);

    case 3

        % Nodes and dofs on the left face of the domain
        % (left and right are defined along the x-direction)
        Mesh.left_nodes = find(Mesh.coords(:,1)==min(Mesh.coords(:,1)));
        Mesh.left_dof = Mesh.DOF(Mesh.left_nodes,:,:);
        Mesh.left_dofx = Mesh.left_dof(:,1);
        Mesh.left_dofy = Mesh.left_dof(:,2);
        Mesh.left_dofz = Mesh.left_dof(:,3);
        Mesh.left_dof = reshape(Mesh.left_dof',numel(Mesh.left_dof),1);

        % Nodes and dofs on the right face of the domain
        Mesh.right_nodes = find(Mesh.coords(:,1)==max(Mesh.coords(:,1)));
        Mesh.right_dof = Mesh.DOF(Mesh.right_nodes,:);
        Mesh.right_dofx = Mesh.right_dof(:,1);
        Mesh.right_dofy = Mesh.right_dof(:,2);
        Mesh.right_dofz = Mesh.right_dof(:,3);
        Mesh.right_dof = reshape(Mesh.right_dof',numel(Mesh.right_dof),1);

        % Nodes and dofs on the top face of the domain
        % (top and bottom are defined along the z-direction)
        Mesh.top_nodes = find(Mesh.coords(:,2)==max(Mesh.coords(:,2)));
        Mesh.top_dof = Mesh.DOF(Mesh.top_nodes,:);
        Mesh.top_dofx = Mesh.top_dof(:,1);
        Mesh.top_dofy = Mesh.top_dof(:,2);
        Mesh.top_dofz = Mesh.top_dof(:,3);
        Mesh.top_dof = reshape(Mesh.top_dof',numel(Mesh.top_dof),1);

        % Nodes and dofs on the bottom face of the domain
        Mesh.bottom_nodes = find(Mesh.coords(:,2)==min(Mesh.coords(:,2)));
        Mesh.bottom_dof =  Mesh.DOF(Mesh.bottom_nodes,:);
        Mesh.bottom_dofx = Mesh.bottom_dof(:,1);
        Mesh.bottom_dofy = Mesh.bottom_dof(:,2);
        Mesh.bottom_dofz = Mesh.bottom_dof(:,3);
        Mesh.bottom_dof = reshape(Mesh.bottom_dof',numel(Mesh.bottom_dof),1);

        % Nodes and dofs on the near face of the domain
        % (near and far are defined along the y-direction)
        Mesh.near_nodes = find(Mesh.coords(:,3)==min(Mesh.coords(:,3)));
        Mesh.near_dof =  Mesh.DOF(Mesh.near_nodes,:);
        Mesh.near_dofx = Mesh.near_dof(:,1);
        Mesh.near_dofy = Mesh.near_dof(:,2);
        Mesh.near_dofz = Mesh.near_dof(:,3);
        Mesh.near_dof = reshape(Mesh.near_dof',numel(Mesh.near_dof),1);

        % Nodes and dofs on the far face of the domain
        Mesh.far_nodes = find(Mesh.coords(:,3)==max(Mesh.coords(:,3)));
        Mesh.far_dof =  Mesh.DOF(Mesh.far_nodes,:);
        Mesh.far_dofx = Mesh.far_dof(:,1);
        Mesh.far_dofy = Mesh.far_dof(:,2);
        Mesh.far_dofz = Mesh.far_dof(:,3);
        Mesh.far_dof = reshape(Mesh.far_dof',numel(Mesh.far_dof),1);

        Mesh.xdofs = 1:3:Mesh.nDOF;
        Mesh.ydofs = 2:3:Mesh.nDOF;
        Mesh.zdofs = 3:3:Mesh.nDOF;
end


end