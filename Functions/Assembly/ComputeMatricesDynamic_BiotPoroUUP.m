function [Kss, Ksp, Mss, Csf, Css, Kpf, Kps, Kpp, Kfp, Mff, Cff, Cfs] = ComputeMatricesDynamic_BiotPoroUUP(Material, MeshU, MeshP, QuadU, QuadP)
% ------------------------------------------------------------------------
% Compute System Matrices for dynamic simulation
% ------------------------------------------------------------------------
% Input parameters: Material, Mesh, Control, Quad
% Output matrices: Kss, Ksp, Ksf, Kpf, Kps, Kpp, Kfs, Kfp, Kff, Mss, Mff
% ------------------------------------------------------------------------

ne = MeshU.ne; % number of elements
nqU = QuadU.nq; % total number of integration points
nqP = QuadP.nq;

% constitutive matrix
C = getConstitutiveMatrix(Material, MeshU);

%% Initialize global matrices
% initialize vector sizes
rowu = zeros(ne*MeshU.nDOFe^2,1);
colu = zeros(ne*MeshU.nDOFe^2,1);

rowp = zeros(ne*MeshP.nDOFe^2,1);
colp = zeros(ne*MeshP.nDOFe^2,1);

rowup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
colup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

Mssvec = zeros(ne*MeshU.nDOFe^2,1);
Mffvec = zeros(ne*MeshU.nDOFe^2,1);

Cssvec = zeros(ne*MeshU.nDOFe^2,1);
Csfvec = zeros(ne*MeshU.nDOFe*MeshU.nDOFe,1);
Cfsvec = zeros(ne*MeshU.nDOFe*MeshU.nDOFe,1);
Cffvec = zeros(ne*MeshU.nDOFe^2,1);

Kssvec = zeros(ne*MeshU.nDOFe^2,1);
Kspvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
Kppvec = zeros(ne*MeshP.nDOFe^2,1);
Kfpvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

% DOF counter
count_u = 1;
count_p = 1;
count_up = 1;

%% Coupled matrices
for e = 1:ne
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
    Mss_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Mff_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    
    Css_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Csf_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Cfs_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Cff_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    
    Kss_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Ksp_e = zeros(MeshU.nDOFe, MeshP.nDOFe);
    Kfp_e = zeros(MeshU.nDOFe, MeshP.nDOFe);
    Kpp_e = zeros(MeshP.nDOFe, MeshP.nDOFe);

    % loop over integration points - Displacement
    for ip = 1:nqU

        % N matrices
        Nu = getN(MeshU, QuadU, ip);

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
        BuVoigt = getBVoigt(MeshU, Bu);

        % assemble local matrices
        Kss_e = Kss_e + (BuVoigt.') * C * BuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        
        Css_e = Css_e + Material.n^2/Material.kf * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Csf_e = Csf_e + Material.n^2/Material.kf * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1) + ...
            (Material.n*Material.mu + Material.n*Material.xif + Material.n*Material.mu/3 - Material.xif*Material.deltaF) * ...
            (BuVoigt.') * BuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Cff_e = Cff_e + Material.n^2/Material.kf * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1) + ...
            (Material.n*Material.mu + Material.n*Material.xif + Material.n*Material.mu/3 - Material.xif*Material.deltaF) * ...
            (BuVoigt.') * BuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Cfs_e = Cfs_e + Material.n^2/Material.kf * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1) - ...
            Material.xif*Material.deltaS * (BuVoigt.') * BuVoigt * Material.t * Jdet * QuadU.w(ip,1);

        Mss_e = Mss_e + (1-Material.n) * Material.rho_s * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Mff_e = Mff_e + Material.n * Material.rho_f * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1);
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

        % changing matrices to Voigt form
        NpVoigt = getNVoigt(MeshP, Np);
        BuVoigt = getBVoigt(MeshU, Bu);

        % assemble local matrices
        Kpp_e = Kpp_e + Material.Minv * (NpVoigt.') * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);

        if MeshU.nsd == 2
            m = [1; 1; 0]; % mapping vector for plane stress
            Ksp_e = Ksp_e + (Material.alpha - Material.n) * (BuVoigt.') * m * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);
            Kfp_e = Kfp_e + Material.n * (BuVoigt.') * m * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);
        else
            Ksp_e = Ksp_e + (Material.alpha - Material.n) * (BuVoigt.') * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);
            Kfp_e = Kfp_e + Material.n * (BuVoigt.') * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);
        end
    end

    % lumped element mass matrix
    if isfield(Material, 'lumpedMass')
        Mss_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        Mff_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        for k = 1:MeshU.nDOFe
            Mss_eDiag(k,k) = sum(Mss_e(k,:));
            Mff_eDiag(k,k) = sum(Mff_e(k,:));
        end
        Mss_e = Mss_eDiag;
        Mff_e = Mff_eDiag;
    end

    % lumped element damping matrix
    if isfield(Material, 'lumpedDamping')
        Css_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        Cff_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        Csf_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        Cfs_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        for k = 1:MeshU.nDOFe
            Css_eDiag(k,k) = sum(Css_e(k,:));
            Cff_eDiag(k,k) = sum(Cff_e(k,:));
            Csf_eDiag(k,k) = sum(Csf_e(k,:));
            Cfs_eDiag(k,k) = sum(Cfs_e(k,:));
        end
        Css_e = Css_eDiag;
        Cff_e = Cff_eDiag;
        Csf_e = Csf_eDiag;
        Cfs_e = Cfs_eDiag;
    end
    
    % vectorized matrices
    count_u = count_u + MeshU.nDOFe^2;
    count_p = count_p + MeshP.nDOFe^2;
    count_up = count_up + MeshU.nDOFe*MeshP.nDOFe;
    
    Mss_e = reshape(Mss_e, [MeshU.nDOFe^2,1]);
    Mff_e = reshape(Mff_e, [MeshU.nDOFe^2,1]);
    
    Css_e = reshape(Css_e, [MeshU.nDOFe^2,1]);
    Csf_e = reshape(Csf_e, [MeshU.nDOFe^2,1]);
    Cfs_e = reshape(Cfs_e, [MeshU.nDOFe^2,1]);
    Cff_e = reshape(Cff_e, [MeshU.nDOFe^2,1]);
    
    Kss_e = reshape(Kss_e, [MeshU.nDOFe^2,1]);
    Ksp_e = reshape(Ksp_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Kfp_e = reshape(Kfp_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Kpp_e = reshape(Kpp_e, [MeshP.nDOFe^2,1]);
    
    rowmatrix_u = dofu_e*ones(1,MeshU.nDOFe);
    rowu_e = reshape(rowmatrix_u, [MeshU.nDOFe^2,1]);
    colu_e = reshape(rowmatrix_u', [MeshU.nDOFe^2,1]);

    rowmatrix_p = dofp_e*ones(1,MeshP.nDOFe);
    rowp_e = reshape(rowmatrix_p, [MeshP.nDOFe^2,1]);
    colp_e = reshape(rowmatrix_p', [MeshP.nDOFe^2,1]);
    
    rowup_e = reshape(dofu_e*ones(1,MeshP.nDOFe),[MeshU.nDOFe*MeshP.nDOFe,1]);
    colup_e = reshape(ones(MeshU.nDOFe,1)*dofp_e.',[MeshU.nDOFe*MeshP.nDOFe,1]);

    Mssvec(count_u-MeshU.nDOFe^2:count_u-1) = Mss_e;
    Mffvec(count_u-MeshU.nDOFe^2:count_u-1) = Mff_e;

    Cssvec(count_u-MeshU.nDOFe^2:count_u-1) = Css_e;
    Csfvec(count_u-MeshU.nDOFe^2:count_u-1) = Csf_e;
    Cfsvec(count_u-MeshU.nDOFe^2:count_u-1) = Cfs_e;
    Cffvec(count_u-MeshU.nDOFe^2:count_u-1) = Cff_e;

    Kssvec(count_u-MeshU.nDOFe^2:count_u-1) = Kss_e;
    Kspvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Ksp_e;
    Kfpvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Kfp_e;
    Kppvec(count_p-MeshP.nDOFe^2:count_p-1) = Kpp_e;

    rowu(count_u-MeshU.nDOFe^2:count_u-1) = rowu_e;
    colu(count_u-MeshU.nDOFe^2:count_u-1) = colu_e;

    rowp(count_p-MeshP.nDOFe^2:count_p-1) = rowp_e;
    colp(count_p-MeshP.nDOFe^2:count_p-1) = colp_e;

    rowup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = rowup_e;
    colup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = colup_e;
end

% sparse matrices
Mss = sparse(rowu, colu, Mssvec, MeshU.nDOF, MeshU.nDOF);
Mff = sparse(rowu, colu, Mffvec, MeshU.nDOF, MeshU.nDOF);

Css = sparse(rowu, colu, Cssvec, MeshU.nDOF, MeshU.nDOF);
Csf = sparse(rowu, colu, Csfvec, MeshU.nDOF, MeshU.nDOF);
Cfs = sparse(rowu, colu, Cfsvec, MeshU.nDOF, MeshU.nDOF);
Cff = sparse(rowu, colu, Cffvec, MeshU.nDOF, MeshU.nDOF);

Kss = sparse(rowu, colu, Kssvec, MeshU.nDOF, MeshU.nDOF);
Ksp = sparse(rowup, colup, Kspvec, MeshU.nDOF, MeshP.nDOF);
Kfp = sparse(rowup, colup, Kfpvec, MeshU.nDOF, MeshP.nDOF);
Kpp = sparse(rowp, colp, Kppvec, MeshP.nDOF, MeshP.nDOF);
Kps = Ksp.';
Kpf = Kfp.';

end