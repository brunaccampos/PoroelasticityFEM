function Mesh = Build1DMesh(nsd, ne, L, type)
% Build Mesh function
% Input parameters: ne = number of elements
%                   L = mesh size
%                   type = element type: L2, L3
% Assumption of uniform mesh

%% Mesh properties
Mesh.ne = ne; % number of elements
Mesh.L = L; % column length
Mesh.nsd = nsd; % number of spatial dimensions
Mesh.type = type;

switch type
    % pressure, porosity
    case 'L2'
        % number of nodes per element
        Mesh.nne = 2;
        % number of DOFs for p per element
        Mesh.nDOFe = 2;
        % number of DOFs for p
        Mesh.nDOF = ne * (Mesh.nDOFe - 1) + 1;
        % DOFs
        Mesh.DOF = (1:Mesh.nDOF).';
        % nodal coordinates for pressure field
        Mesh.coords = (0:(L/(Mesh.nDOF-1)):L).';
        % total number of nodes
        Mesh.nn = length(Mesh.coords);
        % connectivity for dpressure field
        Mesh.conn = zeros(Mesh.ne, Mesh.nDOFe);
        for e = 1:Mesh.ne
            if e == 1
                Mesh.conn(1,:) = 1:Mesh.nDOFe;
            else
                Mesh.conn(e,:) = Mesh.conn(e-1,end):(Mesh.conn(e-1,end) + Mesh.nDOFe -1);
            end
        end

    % displacement
    case 'L3'
        % number of nodes per element
        Mesh.nne = 3;
        % number of DOFs for u per element
        Mesh.nDOFe = 3;
        % number of DOFs for u
        Mesh.nDOF = ne * (Mesh.nDOFe - 1) + 1;
        % DOFs
        Mesh.DOF = (1:Mesh.nDOF).';
        % nodal coordinates for displacement field
        Mesh.coords = (0:(L/(Mesh.nDOF-1)):L).';
        % total number of nodes
        Mesh.nn = length(Mesh.coords);
        % connectivity for displacements field
        Mesh.conn = zeros(Mesh.ne, Mesh.nDOFe);
        for e = 1:Mesh.ne
            if e == 1
                Mesh.conn(1,:) = 1:Mesh.nDOFe;
            else
                Mesh.conn(e,:) = Mesh.conn(e-1,end):(Mesh.conn(e-1,end) + Mesh.nDOFe -1);
            end
        end

end

end