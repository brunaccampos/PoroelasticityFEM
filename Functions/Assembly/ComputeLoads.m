function [fu,fp,fn] = ComputeLoads(BC, MeshU, MeshP, MeshN, Control, QuadU, QuadP, Material)
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

% ------------------------------------------------------------------------
% version 2: different number of IPs for u-p fields
% ------------------------------------------------------------------------

ne = MeshU.ne; % number of elements
nqU = QuadU.nq; % total number of integration points
nqP = QuadP.nq;

% global vectors
fu = zeros(MeshU.nDOF, 1);
fp = zeros(MeshP.nDOF, 1);
if contains(Control.Biotmodel, 'Spanos')
    fn = zeros(MeshN.nDOF, 1);
else
    fn = [];
end

%% Traction vector
% veryfing if there are applied loads
if ~isempty(BC.tractionNodes)
    % traction nodes
    tractionNodes = BC.tractionNodes;
    % traction values
    tractionForce = BC.tractionForce;
    % counter
    count = zeros(size(tractionNodes));
else
    tractionNodes = [];
end

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
    fb_e = zeros(MeshU.nDOFe,1);

    % loop over IPs if there are body forces
    if ~strcmp(func2str(BC.b),'@(x)[]')
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
            
            % body force
            fb_e = fb_e + NVoigt.' * BC.b(ipcoords) * Jdet * QuadU.w(ip,1);
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
            else
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

%     if Control.t <= 1
%         fu_e = fu_e*Control.t;
%     end
    % assemble global load vector
    fu(dofe) = fu(dofe) + fu_e + fb_e;
end

% if Control.t <= 1
%     BC.pointLoad = BC.pointLoad*Control.t;
% end
    
% BC.pointLoad = BC.pointLoadValue * (1 - cos(75*Control.t));

% adding point loads
if ~isempty(BC.pointLoad)
    fu = fu + BC.pointLoad;
end

%% Flux vector
% veryfing if there are applied fluxed
if ~isempty(BC.fluxNodes)
    % traction nodes
    fluxNodes = BC.fluxNodes;
    % traction values
    fluxValue = BC.fluxValue;
    % counter
    count = zeros(size(fluxNodes));
else
    fluxNodes = [];
end

% loop over elements
for e = 1:ne
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
    fs_e = zeros(MeshP.nDOFe,1);

    % loop over IPs if there are flux sources
    if ~strcmp(func2str(BC.s),'@(x)[]')
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
            fs_e = fs_e + NVoigt.' * BC.s(ipcoords) * Jdet * QuadP.w(ip,1);
        end
    end
    
    % loop over flux values
    for j = 1:length(fluxNodes)
        % check if flux is applied to current element
        if ~isempty(find(connp_e == fluxNodes(j),1)) && count(j) == 0
            % current node
            node = find(connp_e == fluxNodes(j),1);
            % node parent coordinates
            coord = getParentCoords(node, MeshP.type);
            % node shape functions
            [N,~] = lagrange_basis(MeshP, coord);
            Nvoigt = getNVoigt(MeshP, N');
            % element load vector
            fp_e = fp_e - Nvoigt.' * fluxValue(j,:)';
            % update counting
            count(j) = 1;
        end
    end

    % assemble global flux vector
    fp(dofe) = fp(dofe) + fp_e + fs_e;
end

% if Control.t <= 100*1e-4
%     BC.pointFlux = BC.pointFlux*Control.t/(100*1e-4);
% end

% adding point loads
if ~isempty(BC.pointFlux)
    fp = fp - BC.pointFlux;
end
    
end