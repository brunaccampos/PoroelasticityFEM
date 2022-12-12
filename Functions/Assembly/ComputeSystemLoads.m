function [fu,fp] = ComputeSystemLoads(BC, MeshU, MeshP, Control, Quad)
% Compute system load force vectors
% ------------------------------------------------------------------------
%   Input
% ------------------------------------------------------------------------
%   BC:     structure array with fields
%       BC.tractionNodes: nodes where traction is applied
%       BC.tractionForce: magnitude of applied traction
%       BC.pointLoad: point loads applied
% ------------------------------------------------------------------------
%   Output
% ------------------------------------------------------------------------
%   fu: global load vector for displacement field
%   fp: global load vector for pressure field (zero vector so far)
% ------------------------------------------------------------------------
% Adapted from: https://github.com/GCMLab (Acknowledgements: Chris Ladubec)
% ------------------------------------------------------------------------

ne = MeshU.ne; % number of elements
nq = Control.nq^MeshU.nsd; % total number of integration points

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

%% Initialize global matrices
fu = zeros(MeshU.nDOF, 1);
fp = zeros(MeshP.nDOF, 1);

%% Load force vector displacement field
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
        for ip = 1:nq
            % N matrices
            N = getN(MeshU, Quad, ip);
            % N derivatives
            dN = getdN(MeshU, Quad, ip);
            % change to Voigt form
            NVoigt = getNVoigt(MeshU, N);
            % quadrature point in phyisical coordinates
            ipcoords = gcoords' * N.';
            % Jacobian matrix
            J = dN*gcoords;
            % Jacobian determinant
            Jdet = det(J);
            
            % body force
            fb_e = fb_e + NVoigt.' * BC.b(ipcoords) * Jdet * Quad.w(ip,1);
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
            fu_e = fu_e + Nvoigt.' * tractionForce(j,:)';
            % update counting
            count(j) = 1;
        end
    end

    % assemble global load vector
    fu(dofe) = fu(dofe) + fu_e + fb_e;
end

% adding point loads
if ~isempty(BC.pointLoad)
    fu = fu + BC.pointLoad;
end