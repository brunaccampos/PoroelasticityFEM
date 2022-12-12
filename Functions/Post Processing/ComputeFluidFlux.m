function [q] = ComputeFluidFlux(Material, Mesh, p)
% Compute strain and stress in the elements of displacement field
% ------------------------------------------------------------------------
% Reference: https://github.com/GCMLab (Acknowledgements: Matin Parchei
% Esfahani)
% ------------------------------------------------------------------------

%% Initialize variables
nn = Mesh.nn; % total number of nodes
ne = Mesh.ne; % number of elements
nne = Mesh.nne; % number of nodes per element
kf = Material.kf; % porous media permeability [m2/Pa s]

% dimension of flux matrix
switch Mesh.nsd
    case 1
        dim = 1;
    case 2
        dim = 2;
end

q = zeros(nn,dim); % flux matrix
count = zeros(nn,dim);

%% Loop over elements
for e = 1:ne
    % element connectivity
    conne = Mesh.conn(e,:);
    % element DOF numbers
    dofe = Mesh.DOF(conne,:);
    dofe = reshape(dofe', Mesh.nDOFe,[]);
    % element global coordinates
    gcoords = Mesh.coords(conne,:);
    % element nodal pressure
    pe = p(dofe);
    
    % element matrices
    q_e = zeros(nne, dim);
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
        % element flux
        q_e(n,:) = -kf * (B*pe).';
    end
    
    % add to global matrices
    q(conne,:) = q(conne,:) + q_e;
    count(conne,:) = count(conne,:) + 1;
end

% average over nodes
q = q./count;

end