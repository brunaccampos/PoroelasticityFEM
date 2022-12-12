function [fu,fp] = ComputeSystemLoads(BC, MeshU, MeshP)
% Compute system load force vectors
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   BC:     structure array with fields
%       BC.tractionNodes: nodes where traction is applied
%       BC.tractionForce: magnitude of applied traction
%       BC.pointLoad: point loads applied
%   --------------------------------------------------------------------
%   Output
%   --------------------------------------------------------------------
%   fu: global load vector for displacement field
%   fp: global load vector for pressure field (zero vector so far)
%
% Adapted from: https://github.com/GCMLab (Acknowledgements: Chris Ladubec)

ne = MeshU.ne; % number of elements

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
    % element DOFs
    dofe = MeshU.DOF(connu_e,:);
    dofe = reshape(dofe', MeshU.nDOFe,[]);
    % element load matrix
    fu_e = zeros(MeshU.nDOFe,1);

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
        end
    end

    % assemble global load vector
    fu(dofe) = fu_e;
end

% adding point loads
if ~isempty(BC.pointLoad)
    fu = fu + BC.pointLoad;
end