function [strainU, stressU] = ComputeSolidStress(Material, Mesh, Control, Quad, u)
% Compute strain and stress in the elements of displacement field

% Reference: https://github.com/GCMLab (Acknowledgements: Matin Parchei
% Esfahani)

%% Initialize variables
nn = Mesh.nn; % total number of nodes
ne = Mesh.ne; % number of elements
nne = Mesh.nne; % number of nodes per element
% nq = Control.nq^Mesh.nsd; % total number of integration points
C = getConstitutiveMatrix(Material, Mesh); % constitutive matrix

strainU = zeros(nn,3); % strain matrix
stressU = zeros(nn,3); % stress matrix
count = zeros(nn,1);

%% Loop over elements
for e = 1:ne
    % element connectivity
    conne = Mesh.conn(e,:);
    % element DOF numbers
    dofe = Mesh.DOF(conne,:);
    dofe = reshape(dofe', Mesh.nDOFe,[]);
    % element global coordinates
    gcoords = Mesh.coords(conne,:);
    % element nodal displacements
    ue = u(dofe);
    
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
        strainU(conne(n),:) = strainU(conne(n),:) + (B*ue).';
        % stress
        stressU(conne(n),:) = stressU(conne(n),:) + (C*strainU(conne(n),:).').';
        
        count(conne(n)) = count(conne(n)) + 1;
    end
end

% average over nodes
strainU(:,1) = strainU(:,1)./count;
strainU(:,2) = strainU(:,2)./count;
strainU(:,3) = strainU(:,3)./count;

stressU(:,1) = stressU(:,1)./count;
stressU(:,2) = stressU(:,2)./count;
stressU(:,3) = stressU(:,3)./count;

end