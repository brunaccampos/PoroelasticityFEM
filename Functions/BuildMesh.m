function Mesh = BuildMesh(ne, L)
% Build Mesh function
% Input parameters: ne = number of elements
%                   L = mesh size
% Assumption of uniform mesh

%% Mesh properties
Mesh.ne = ne;
Mesh.L = L;

%% Degrees of freedom
% number of DOFs for u per element
Mesh.ndof_u_e = 3; 
% number of DOFs for p per element
Mesh.ndof_p_e = 2; 
% number of DOFs for u
Mesh.ndof_u = ne * (Mesh.ndof_u_e - 1) + 1; 
% number of DOFs for p
Mesh.ndof_p = ne * (Mesh.ndof_p_e - 1) + 1; 

% DOFs
Mesh.dof_u = (1:Mesh.ndof_u);
Mesh.dof_p = (1:Mesh.ndof_p);

% element size [m]
Mesh.h = L/ne;

%% Nodal coordinates
% nodal coordinates for displacement field
Mesh.coords_u = (0:(L/(Mesh.ndof_u-1)):L);
% nodal coordinates for pressure field
Mesh.coords_p = (0:(L/(Mesh.ndof_p-1)):L); 

%% Connectivity matrices
% connectivity for displacements field
Mesh.connu = zeros(Mesh.ne, Mesh.ndof_u_e);
for e = 1:Mesh.ne
    if e == 1
        Mesh.connu(1,:) = 1:Mesh.ndof_u_e;
    else
        Mesh.connu(e,:) = Mesh.connu(e-1,end):(Mesh.connu(e-1,end) + Mesh.ndof_u_e -1);
    end
end

% connectivity for dpressure field
Mesh.connp = zeros(Mesh.ne, Mesh.ndof_p_e);
for e = 1:Mesh.ne
    if e == 1
        Mesh.connp(1,:) = 1:Mesh.ndof_p_e;
    else
        Mesh.connp(e,:) = Mesh.connp(e-1,end):(Mesh.connp(e-1,end) + Mesh.ndof_p_e -1);
    end
end

end