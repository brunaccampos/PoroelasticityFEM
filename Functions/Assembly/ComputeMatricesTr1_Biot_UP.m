function [Kuu, Kup, Kpp, Kpu, S] = ComputeMatricesTr1_Biot_UP(Material, MeshU, MeshP, QuadU, QuadP)
% Compute System Matrices for 1D quasi-steady simulation
% Input parameters: Material, Mesh, Control, Quad
% Output matrices: Kuu, Kup, Kpp, S
% ------------------------------------------------------------------------
% version 4: different number of IPs for u-p fields
% ------------------------------------------------------------------------

ne = MeshU.ne; % number of elements
nqU = QuadU.nq; % total number of integration points
nqP = QuadP.nq;

%% Initialize global matrices
% initialize vector sizes
rowu = zeros(ne*MeshU.nDOFe^2,1);
colu = zeros(ne*MeshU.nDOFe^2,1);

rowp = zeros(ne*MeshP.nDOFe^2,1);
colp = zeros(ne*MeshP.nDOFe^2,1);

rowup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
colup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

Kuuvec = zeros(ne*MeshU.nDOFe^2,1);
Kppvec = zeros(ne*MeshP.nDOFe^2,1);
Svec = zeros(ne*MeshP.nDOFe^2,1);
Kupvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

% DOF counter
count_u = 1;
count_p = 1;
count_up = 1;

%% Coupled matrices
for e = 1:ne
    % element material type
    nMat = MeshU.MatList(e); % element material type
    % constitutive matrix
    C = getConstitutiveMatrix(nMat, Material, MeshU);

    % element connectivity
    connu_e = MeshU.conn(e,:);
    connu_e = reshape(connu_e',MeshU.nne,[]);
    connp_e = MeshP.conn(e,:);
    connp_e = reshape(connp_e',MeshP.nne,[]);

    % element DOF numbers
    dofu_e = MeshU.DOF(connu_e,:);
    dofu_e = reshape(dofu_e',MeshU.nDOFe,[]);
    dofp_e = MeshP.DOF(connp_e,:);
    dofp_e = reshape(dofp_e',MeshP.nDOFe,[]);

    % global coordinates
    gcoordsU = MeshU.coords(connu_e,:);
    gcoordsP = MeshP.coords(connp_e,:);

    % initialize local matrices
    Kuu_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Kup_e = zeros(MeshU.nDOFe,MeshP.nDOFe);
    Kpp_e = zeros (MeshP.nDOFe, MeshP.nDOFe);
    S_e = zeros (MeshP.nDOFe, MeshP.nDOFe);

    % loop over integration points
    for ip = 1:nqU

        % N matrices
        Np = getN(MeshP, QuadU, ip);

        % N derivatives
        dNu = getdN(MeshU, QuadU, ip);
        dNp = getdN(MeshP, QuadU, ip);

        % Jacobian matrix
        Ju = dNu*gcoordsU;
        Jp = dNp*gcoordsP;
        % Jacobian determinant
        Jdet = det(Jp);

        % B matrices
        Bu = Ju\dNu;

        % changing matrices to Voigt form
        NpVoigt = getNVoigt(MeshP, Np);
        BuVoigt = getBVoigt(MeshU, Bu);

        % assemble local matrices
        Kuu_e = Kuu_e + (BuVoigt.') * C * BuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        S_e = S_e + Material.M(nMat).Minv * (NpVoigt.') * NpVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Kup_e = Kup_e + Material.M(nMat).alpha * (BuVoigt.') * Material.m * NpVoigt * Material.t * Jdet * QuadU.w(ip,1);
    end
    
    % loop over integration points
    for ip = 1:nqP

        % N derivatives
        dNp = getdN(MeshP, QuadP, ip);

        % Jacobian matrix
        Jp = dNp*gcoordsP;
        % Jacobian determinant
        Jdet = det(Jp);

        % B matrices
        Bp = Jp\dNp;

        % changing matrices to Voigt form
        BpVoigt = getBVoigt(MeshP, Bp);

        % assemble local matrices
        Kpp_e = Kpp_e + Material.M(nMat).kf * (BpVoigt.') * BpVoigt * Material.t * Jdet * QuadP.w(ip,1);        
    end
    
    % vectorized matrices
    count_u = count_u + MeshU.nDOFe^2;
    count_p = count_p + MeshP.nDOFe^2;
    count_up = count_up + MeshU.nDOFe*MeshP.nDOFe;
    
    Kuu_e = reshape(Kuu_e, [MeshU.nDOFe^2,1]);
    Kpp_e = reshape(Kpp_e, [MeshP.nDOFe^2,1]);
    S_e = reshape(S_e, [MeshP.nDOFe^2,1]);
    Kup_e = reshape(Kup_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    
    rowmatrix_u = dofu_e*ones(1,MeshU.nDOFe);
    rowu_e = reshape(rowmatrix_u, [MeshU.nDOFe^2,1]);
    colu_e = reshape(rowmatrix_u', [MeshU.nDOFe^2,1]);

    rowmatrix_p = dofp_e*ones(1,MeshP.nDOFe);
    rowp_e = reshape(rowmatrix_p, [MeshP.nDOFe^2,1]);
    colp_e = reshape(rowmatrix_p', [MeshP.nDOFe^2,1]);

    rowup_e = reshape(dofu_e*ones(1,MeshP.nDOFe),[MeshU.nDOFe*MeshP.nDOFe,1]);
    colup_e = reshape(ones(MeshU.nDOFe,1)*dofp_e.',[MeshU.nDOFe*MeshP.nDOFe,1]);
    
    Kuuvec(count_u-MeshU.nDOFe^2:count_u-1) = Kuu_e;
    Kppvec(count_p-MeshP.nDOFe^2:count_p-1) = Kpp_e;
    Svec(count_p-MeshP.nDOFe^2:count_p-1) = S_e;
    Kupvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Kup_e;
    
    rowu(count_u-MeshU.nDOFe^2:count_u-1) = rowu_e;
    colu(count_u-MeshU.nDOFe^2:count_u-1) = colu_e;

    rowp(count_p-MeshP.nDOFe^2:count_p-1) = rowp_e;
    colp(count_p-MeshP.nDOFe^2:count_p-1) = colp_e;

    rowup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = rowup_e;
    colup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = colup_e;
end

% sparse matrices
Kuu = sparse(rowu, colu, Kuuvec, MeshU.nDOF, MeshU.nDOF);
Kpp = sparse(rowp, colp, Kppvec, MeshP.nDOF, MeshP.nDOF);
S = sparse(rowp, colp, Svec, MeshP.nDOF, MeshP.nDOF);
Kup = sparse(rowup, colup, Kupvec, MeshU.nDOF, MeshP.nDOF);
Kpu = Kup.';

end