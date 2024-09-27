function [Kuu, Kup, Kpp, Muu, Mpu, S] = ComputeMatricesDyn_BT_UP(Material, MeshU, MeshP, QuadU, QuadP)
% ------------------------------------------------------------------------
% Compute System Matrices for 1D dynamic simulation
% ------------------------------------------------------------------------
% Input parameters: Material, Mesh, Control, Quad
% Output matrices: Kuu, Kup, Kpp, M, Mhat, S
% ------------------------------------------------------------------------

ne = MeshU.ne; % number of elements
nqU = QuadU.nq; % total number of integration points
nqP = QuadP.nq;

%% Initialize global matrices
% initialize vector sizes
% u-u
rowu = zeros(ne*MeshU.nDOFe^2,1);
colu = zeros(ne*MeshU.nDOFe^2,1);
% p-p
rowp = zeros(ne*MeshP.nDOFe^2,1);
colp = zeros(ne*MeshP.nDOFe^2,1);
% u-p
rowup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
colup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
% p-u
rowpu = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
colpu = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

% initialize vectorized global matrices
Muuvec = zeros(ne*MeshU.nDOFe^2,1);
Mpuvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

Kuuvec = zeros(ne*MeshU.nDOFe^2,1);
Kppvec = zeros(ne*MeshP.nDOFe^2,1);
Kupvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

Svec = zeros(ne*MeshP.nDOFe^2,1);

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
    Muu_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Mpu_e = zeros(MeshP.nDOFe, MeshU.nDOFe);
    
    Kuu_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Kup_e = zeros(MeshU.nDOFe, MeshP.nDOFe);
    Kpp_e = zeros(MeshP.nDOFe, MeshP.nDOFe);
    
    S_e = zeros(MeshP.nDOFe, MeshP.nDOFe);

    % loop over integration points - Displacement
    for ip = 1:nqU

        % N matrices
        Nu = getN(MeshU, QuadU, ip);
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
        NuVoigt = getNVoigt(MeshU, Nu);
        NpVoigt = getNVoigt(MeshP, Np);
        BuVoigt = getBVoigt(MeshU, Bu);

        % assemble local matrices
        Kuu_e = Kuu_e + (BuVoigt.') * C * BuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Muu_e = Muu_e + Material.M(nMat).rho * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        S_e = S_e + Material.M(nMat).Minv * (NpVoigt.') * NpVoigt * Material.t * Jdet * QuadU.w(ip,1);
    end
 
    % loop over integration points - Displacement
    for ip = 1:nqP

        % N matrices
        Np = getN(MeshP, QuadP, ip);

        % N derivatives
        dNu = getdN(MeshU, QuadP, ip);
        dNp = getdN(MeshP, QuadP, ip);

        % Jacobian matrix
        Ju = dNu*gcoordsU;
        Jp = dNp*gcoordsP;
        % Jacobian determinant
        Jdet = det(Jp);

        % B matrices
        Bu = Ju\dNu;
        Bp = Jp\dNp;

        % changing matrices to Voigt form
        NpVoigt = getNVoigt(MeshP, Np);
        BuVoigt = getBVoigt(MeshU, Bu);
        BpVoigt = getBVoigt(MeshP, Bp);


        % assemble local matrices
        Kpp_e = Kpp_e + Material.M(nMat).kf * (BpVoigt.') * BpVoigt * Material.t * Jdet * QuadP.w(ip,1);
        Kup_e = Kup_e + Material.M(nMat).alpha * (BuVoigt.') * Material.m * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);
        Mpu_e = Mpu_e + Material.M(nMat).rhof * Material.M(nMat).kf * (NpVoigt.') * (Material.m') * BuVoigt * Material.t * Jdet * QuadP.w(ip,1);
    end

    % lumped element mass matrix
    if Material.lumpedMass
        Muu_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        for k = 1:MeshU.nDOFe
            Muu_eDiag(k,k) = sum(Muu_e(k,:));
        end
        Muu_e = Muu_eDiag;
    end
    
    % vectorized matrices
    count_u = count_u + MeshU.nDOFe^2;
    count_p = count_p + MeshP.nDOFe^2;
    count_up = count_up + MeshU.nDOFe*MeshP.nDOFe;
    
    Muu_e = reshape(Muu_e, [MeshU.nDOFe^2,1]);
    Mpu_e = reshape(Mpu_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    
    Kuu_e = reshape(Kuu_e, [MeshU.nDOFe^2,1]);
    Kpp_e = reshape(Kpp_e, [MeshP.nDOFe^2,1]);
    Kup_e = reshape(Kup_e, [MeshU.nDOFe*MeshP.nDOFe,1]);

    S_e = reshape(S_e, [MeshP.nDOFe^2,1]);
    
    % u-u
    rowmatrix_u = dofu_e*ones(1,MeshU.nDOFe);
    rowu_e = reshape(rowmatrix_u, [MeshU.nDOFe^2,1]);
    colu_e = reshape(rowmatrix_u', [MeshU.nDOFe^2,1]);
    % p-p
    rowmatrix_p = dofp_e*ones(1,MeshP.nDOFe);
    rowp_e = reshape(rowmatrix_p, [MeshP.nDOFe^2,1]);
    colp_e = reshape(rowmatrix_p', [MeshP.nDOFe^2,1]);
    % u-p
    rowup_e = reshape(dofu_e*ones(1,MeshP.nDOFe),[MeshU.nDOFe*MeshP.nDOFe,1]);
    colup_e = reshape(ones(MeshU.nDOFe,1)*dofp_e.',[MeshU.nDOFe*MeshP.nDOFe,1]);
    % p-u
    rowpu_e = reshape(dofp_e*ones(1,MeshU.nDOFe),[MeshU.nDOFe*MeshP.nDOFe,1]);
    colpu_e = reshape(ones(MeshP.nDOFe,1)*dofu_e.',[MeshU.nDOFe*MeshP.nDOFe,1]);
    
    Muuvec(count_u-MeshU.nDOFe^2:count_u-1) = Muu_e;
    Mpuvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Mpu_e;
    
    Kuuvec(count_u-MeshU.nDOFe^2:count_u-1) = Kuu_e;
    Kppvec(count_p-MeshP.nDOFe^2:count_p-1) = Kpp_e;
    Kupvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Kup_e;
    
    Svec(count_p-MeshP.nDOFe^2:count_p-1) = S_e;
    
    % u-u
    rowu(count_u-MeshU.nDOFe^2:count_u-1) = rowu_e;
    colu(count_u-MeshU.nDOFe^2:count_u-1) = colu_e;
    % p-p
    rowp(count_p-MeshP.nDOFe^2:count_p-1) = rowp_e;
    colp(count_p-MeshP.nDOFe^2:count_p-1) = colp_e;
    % u-p
    rowup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = rowup_e;
    colup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = colup_e;
    % p-u
    rowpu(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = rowpu_e;
    colpu(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = colpu_e;
end

% sparse matrices
Muu = sparse(rowu, colu, Muuvec, MeshU.nDOF, MeshU.nDOF);
Mpu = sparse(rowpu, colpu, Mpuvec, MeshP.nDOF, MeshU.nDOF);

Kuu = sparse(rowu, colu, Kuuvec, MeshU.nDOF, MeshU.nDOF);
Kpp = sparse(rowp, colp, Kppvec, MeshP.nDOF, MeshP.nDOF);
Kup = sparse(rowup, colup, Kupvec, MeshU.nDOF, MeshP.nDOF);

S = sparse(rowp, colp, Svec, MeshP.nDOF, MeshP.nDOF);

end