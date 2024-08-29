function [fs, fp, ff] = ComputeLoads_UPU(BC, MeshU, MeshP, Control, Material, QuadU, QuadP)
% Compute system load force vectors
% ------------------------------------------------------------------------
%   Input
% ------------------------------------------------------------------------
%   BC:     structure array with fields
%       BC.tractionNodes: nodes where traction is applied
%       BC.tractionForce: magnitude of applied traction
%       BC.pointLoad: point loads applied
%       BC.b: body forces
%
%       BC.fluxNodes: nodes where flux is applied
%       BC.fluxValue: magnitude of applied flux
%       BC.pointFlux: point fluxes applied
%       BC.s: flux sources
% ------------------------------------------------------------------------
%   Output
% ------------------------------------------------------------------------
%   fu: global traction vector for displacement field
%   fp: global flux vector for pressure field
% ------------------------------------------------------------------------
% Adapted from: https://github.com/GCMLab (Acknowledgements: Chris Ladubec)
% ------------------------------------------------------------------------

ne = MeshU.ne; % number of elements
nqU = QuadU.nq; % total number of integration points
nqP = QuadP.nq;

% global vectors
fs = zeros(MeshU.nDOF, 1);
fp = zeros(MeshP.nDOF, 1);
ff = zeros(MeshU.nDOF, 1);

%% Return zeros if no applied loads
% verifying if there are body forces
if isempty(BC.tractionNodes) && strcmp(func2str(BC.bs),'@(x,t)[]') && strcmp(func2str(BC.bf),'@(x,t)[]') && ...
    strcmp(func2str(BC.pointLoad),'@(t)[]') && isempty(BC.fluxNodes) && strcmp(func2str(BC.s),'@(x,t)[]') && strcmp(func2str(BC.pointFlux),'@(t)[]')
    return
end

% verifying if there are distributed forces
if ~isempty(BC.tractionNodes)
    % traction nodes
    tractionNodes = BC.tractionNodes;
    % traction values
    tractionForce = BC.tractionForce(Control.t);
    % counter
    count = zeros(size(tractionNodes));
else
    tractionNodes = [];
end

%% Body force for us-uf fields, traction for us field
% loop over elements
for e = 1:ne
    % element connectivity
    connu_e = MeshU.conn(e,:);
    % global coordinates
    gcoords = MeshU.coords(connu_e,:);
    % element DOFs
    dofe = MeshU.DOF(connu_e,:);
    dofe = reshape(dofe', MeshU.nDOFe,[]);
    % element load matrix
    fu_e = zeros(MeshU.nDOFe,1);
    % element body force matrix
    fbs_e = zeros(MeshU.nDOFe,1); % solid eq.
    fbf_e = zeros(MeshU.nDOFe,1); % fluid eq.
    
    % loop over IPs if there are body forces - solid equation
    if ~strcmp(func2str(BC.bs),'@(x,t)[]')
        for ip = 1:nqU
            % N matrices
            N = getN(MeshU, QuadU, ip);
            % N derivatives
            dN = getdN(MeshU, QuadU, ip);
            % change to Voigt form
            NVoigt = getNVoigt(MeshU, N);
            % quadrature point in phyisical coordinates
            ipcoords = gcoords' * N.';
            % Jacobian matrix
            J = dN*gcoords;
            % Jacobian determinant
            Jdet = det(J);
            
            % body force solid equation
            fbs_e = fbs_e + NVoigt.' * BC.bs(ipcoords, Control.t) * Jdet * QuadU.w(ip,1);
        end
    end
    
    % loop over IPs if there are body forces - fluid equation
    if ~strcmp(func2str(BC.bf),'@(x,t)[]')
        for ip = 1:nqU
            % N matrices
            N = getN(MeshU, QuadU, ip);
            % N derivatives
            dN = getdN(MeshU, QuadU, ip);
            % change to Voigt form
            NVoigt = getNVoigt(MeshU, N);
            % quadrature point in phyisical coordinates
            ipcoords = gcoords' * N.';
            % Jacobian matrix
            J = dN*gcoords;
            % Jacobian determinant
            Jdet = det(J);
            
            % body force fluid equation
            fbf_e = fbf_e + NVoigt.' * BC.bf(ipcoords, Control.t) * Jdet * QuadU.w(ip,1);
        end
    end
    
    % loop over traction forces
    for j = 1:length(tractionNodes)
        % check if traction is applied to current element
        if ~isempty(find(connu_e == tractionNodes(j),1)) && count(j) == 0
            % current node
            node = find(connu_e == tractionNodes(j),1);
            % node parent coordinates
            coord = getParentCoords(node, MeshU.type);
            % node shape functions
            [N,~] = lagrange_basis(MeshU, coord);
            Nvoigt = getNVoigt(MeshU, N');
            % element load vector
            if ~BC.tractionInterp
                fu_e = fu_e + Nvoigt.' * tractionForce(j,:)';
            else % reorganize vector when interpolating
                tractionForce_e = tractionForce(connu_e,:);
                tractionForce_eVec = zeros(MeshU.nne*2,1);
                tractionForce_eVec(1:2:end) = tractionForce_e(:,1);
                tractionForce_eVec(2:2:end) = tractionForce_e(:,2);
                fu_e = fu_e + Nvoigt.' * Nvoigt * tractionForce_eVec;
            end
            % update counting
            count(j) = 1;
        end
    end

    % apply step load gradualy
    if Control.rampLoad
        if Control.t <= Control.tlim
            fu_e = fu_e*Control.t/Control.tlim;
        end
    end

    % assemble global load vectors
    fs(dofe) = fs(dofe) + fu_e + fbs_e;
    ff(dofe) = ff(dofe) + fbf_e;
end

% apply step load gradualy
if Control.rampLoad
    if Control.t <= Control.tlim
        aux = Control.t/Control.tlim;
        BC.pointLoad = @(t) BC.pointLoad(t) * aux;
    end
end
    
%% Point load u field
if ~strcmp(func2str(BC.pointLoad),'@(t)[]')
    fs = fs + BC.pointLoad(Control.t);
end

%% Source/sink term, p vector
% counter for pressure BC nodes
count = zeros(size(BC.fixed_p));

% loop over elements
for e = 1:ne
    % element material type
    nMat = Mesh.MatList(e);
    % element connectivity
    connp_e = MeshP.conn(e,:);
    % global coordinates
    gcoords = MeshP.coords(connp_e,:);
    % element DOFs
    dofe = MeshP.DOF(connp_e,:);
    dofe = reshape(dofe', MeshP.nDOFe,[]);
    % element flux matrix
    fp_e = zeros(MeshP.nDOFe,1);
    % element flux source matrix
    fg_e = zeros(MeshP.nDOFe,1);
    % auxiliar matrix for p vector
    faux_e = zeros(MeshP.nDOFe,1);
    
    % loop over IPs if there are flux sources
    if ~strcmp(func2str(BC.s),'@(x,t)[]')
        for ip = 1:nqP
            % N matrices
            N = getN(MeshP, QuadP, ip);
            % N derivatives
            dN = getdN(MeshP, QuadP, ip);
            % change to Voigt form
            NVoigt = getNVoigt(MeshP, N);
            % quadrature point in phyisical coordinates
            ipcoords = gcoords' * N.';
            % Jacobian matrix
            J = dN*gcoords;
            % Jacobian determinant
            Jdet = det(J);
            
            % flux source
            fg_e = fg_e + NVoigt.' * BC.s(ipcoords, Control.t) * Jdet * QuadP.w(ip,1);
        end
    end
    
    % loop over prescribed pressure
    for j = 1:length(BC.fixed_p)
        % check if flux is applied to current element
        if ~isempty(find(connp_e == BC.fixed_p(j),1)) && count(j) == 0
            % current node
            node = find(connp_e == BC.fixed_p(j),1);
            % node parent coordinates
            coord = getParentCoords(node, MeshP.type);
            % node shape functions
            [N,~] = lagrange_basis(MeshP, coord);
            Nvoigt = getNVoigt(MeshP, N');
            % element load vector
            vec = BC.fixed_p_value(Control.t);
            faux_e = faux_e + Nvoigt.' * vec(j,:)';
            % update counting
            count(j) = 1;
        end
    end
    
    % assemble global flux vector
    fp(dofe) = fp(dofe) + fg_e;
    % add pressure vector contributions (valid for P2/P1 meshes)
    fs(dofe*2-1) = fs(dofe*2-1) - (Material.M(nMat).alpha-Material.M(nMat).eta0) * faux_e;
    ff(dofe*2-1) = ff(dofe*2-1) - Material.M(nMat).eta0 * faux_e;
    
end
   
end