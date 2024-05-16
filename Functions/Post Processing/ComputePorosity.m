function [eta, etadot] = ComputePorosity(Material, Mesh, Solution, Control)
% Compute porosity based on the solid and fluid displacements values
% ------------------------------------------------------------------------

%% Initialize variables
nn = Mesh.nn; % total number of nodes
ne = Mesh.ne; % number of elements
nne = Mesh.nne; % number of nodes per element

us = Solution.u; % solid displacement
usdot = Solution.udot; % solid velocity

% fluid displacement and velocity
if contains(Control.PMmodel, 'UPU')
    uf = Solution.uf;
    ufdot = Solution.ufdot;
elseif contains(Control.PMmodel, 'UPV')
    uf = zeros(MeshU.nDOF,1);
    ufdot = Solution.ufdot;
elseif contains(Control.PMmodel, 'UPW')
    uf = zeros(MeshU.nDOF,1);
    ufdot = (Solution.w + Material.eta0*Solution.udot)/Material.eta0;
end

eta = zeros(nn,1); % porosity
etadot = zeros(nn,1); % time-varying porosity
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
    % element nodal displacements (solid)
    use = us(dofe);
    % element nodal displacements (fluid)
    ufe = uf(dofe);
    % element nodal velocities (solid)
    usdote = usdot(dofe);
    % element nodal velocities (fluid)
    ufdote = ufdot(dofe);
    
    % element matrices
    eta_e = zeros(nne, 1);
    etadot_e = zeros(nne, 1);
    
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
        % porosity
        eta_e(n,:) = (B*use).' - (B*ufe).' + Material.eta0;
        % time varying porosity
        etadot_e(n,:) = (B*usdote).' - (B*ufdote).';
    end
    
    % add to global matrices
    eta(conne,:) = eta(conne,:) + eta_e;
    etadot(conne,:) = etadot(conne,:) + etadot_e;
    
    % update counter
    count(conne,:) = count(conne,:) + 1;
end

% average over nodes
eta = eta./count;
etadot = etadot./count;

end