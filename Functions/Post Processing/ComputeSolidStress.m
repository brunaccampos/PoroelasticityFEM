% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [strain, stress] = ComputeSolidStress(Material, Mesh, u)
% Compute strain and stress in the elements of displacement field
% Acknowledgements: Matin Parchei Esfahani

%% Initialize variables
nn = Mesh.nn; % total number of nodes
ne = Mesh.ne; % number of elements
nne = Mesh.nne; % number of nodes per element

% dimension of strain/stress matrix
switch Mesh.nsd
    case 1
        dim = 1;
    case 2
        dim = 3;
end

strain = zeros(nn,dim); % strain matrix
stress = zeros(nn,dim); % stress matrix
count = zeros(nn,dim);

%% Loop over elements
for e = 1:ne
    % element material type
    nMat = Mesh.MatList(e);
    % constitutive matrix
    C = getConstitutiveMatrix(nMat, Material, Mesh);
    
    % element connectivity
    conne = Mesh.conn(e,:);
    % element DOF numbers
    dofe = Mesh.DOF(conne,:);
    dofe = reshape(dofe', Mesh.nDOFe,[]);
    % element global coordinates
    gcoords = Mesh.coords(conne,:);
    % element nodal displacements
    ue = u(dofe);
    
    % element matrices
    strain_e = zeros(nne, dim);
    stress_e = zeros(nne, dim);
    % loop over element nodes
    for n = 1:nne
        % node parent coordinates
        coord = getParentCoords(n, Mesh.type);
        % shape function derivatives
        [~,dN] = lagrange_basis(Mesh, coord);
        % Jacobian matrix
        J = dN.'*gcoords;
        % B matrix
        B = J\dN.';
        % changing to Voigt form
        B = getBVoigt(Mesh,B);
        % strain
        strain_e(n,:) = (B*ue).';
        % stress
        stress_e(n,:) = (C*strain_e(n,:).').';
    end
    
    % add to global matrices
    strain(conne,:) = strain(conne,:) + strain_e;
    stress(conne,:) = stress(conne,:) + stress_e;
    
    % update counter
    count(conne,:) = count(conne,:) + 1;
end

% average over nodes
strain = strain./count;
stress = stress./count;

end