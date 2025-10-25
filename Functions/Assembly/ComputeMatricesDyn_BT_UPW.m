% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Mss, Msf, Kss, Ksp, Mfs, Mff, Kfp, Kff, Cps, Kpf, Cpp, Cfs] = ComputeMatricesDyn_BT_UPW(Material, MeshU, MeshP, QuadU, QuadP)
% Compute system matrices for dynamic simulation

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

Mssvec = zeros(ne*MeshU.nDOFe^2,1);
Mffvec = zeros(ne*MeshU.nDOFe^2,1);
Msfvec = zeros(ne*MeshU.nDOFe^2,1);

Cspvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
Cppvec = zeros(ne*MeshP.nDOFe^2,1);

Kssvec = zeros(ne*MeshU.nDOFe^2,1);
Kspvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
Kffvec = zeros(ne*MeshU.nDOFe^2,1);
Kfpvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);

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
    Mss_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Mff_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Msf_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    
    Csp_e = zeros(MeshU.nDOFe, MeshP.nDOFe);
    Cpp_e = zeros(MeshP.nDOFe, MeshP.nDOFe);
    
    Kss_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Ksp_e = zeros(MeshU.nDOFe, MeshP.nDOFe);
    Kfp_e = zeros(MeshU.nDOFe, MeshP.nDOFe);
    Kff_e = zeros(MeshU.nDOFe, MeshU.nDOFe);

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
        Kff_e = Kff_e + real(Material.M(nMat).F_BT) * Material.M(nMat).muf/Material.M(nMat).k * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadP.w(ip,1);
     
        Mss_e = Mss_e + Material.M(nMat).rho * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Mff_e = Mff_e + Material.M(nMat).rhof/Material.M(nMat).eta0 * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Msf_e = Msf_e + Material.M(nMat).rhof * (NuVoigt.') * NuVoigt * Material.t * Jdet * QuadU.w(ip,1);
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
        Ksp_e = Ksp_e + Material.M(nMat).alpha * (BuVoigt.') * Material.m * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);
        Kfp_e = Kfp_e + (BuVoigt.') * Material.m * NpVoigt * Material.t * Jdet * QuadP.w(ip,1);
        Csp_e = Csp_e + Material.M(nMat).alpha * (BuVoigt.') * Material.m * NpVoigt * Material.t * Jdet * QuadU.w(ip,1);
        Cpp_e = Cpp_e + Material.M(nMat).Minv * (NpVoigt.') * NpVoigt * Material.t * Jdet * QuadU.w(ip,1);
    end

    % lumped element mass matrix
    if Material.lumpedMass
        Mss_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        Mff_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        Msf_eDiag = zeros(MeshU.nDOFe, MeshU.nDOFe);
        for k = 1:MeshU.nDOFe
            Mss_eDiag(k,k) = sum(Mss_e(k,:));
            Mff_eDiag(k,k) = sum(Mff_e(k,:));
            Msf_eDiag(k,k) = sum(Msf_e(k,:));
        end
        Mss_e = Mss_eDiag;
        Mff_e = Mff_eDiag;
        Msf_e = Msf_eDiag;
    end
  
    % vectorized matrices
    count_u = count_u + MeshU.nDOFe^2;
    count_p = count_p + MeshP.nDOFe^2;
    count_up = count_up + MeshU.nDOFe*MeshP.nDOFe;
    
    Mss_e = reshape(Mss_e, [MeshU.nDOFe^2,1]);
    Mff_e = reshape(Mff_e, [MeshU.nDOFe^2,1]);
    Msf_e = reshape(Msf_e, [MeshU.nDOFe^2,1]);
    
    Csp_e = reshape(Csp_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Cpp_e = reshape(Cpp_e, [MeshP.nDOFe^2,1]);
    
    Kss_e = reshape(Kss_e, [MeshU.nDOFe^2,1]);
    Ksp_e = reshape(Ksp_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Kfp_e = reshape(Kfp_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Kff_e = reshape(Kff_e, [MeshU.nDOFe^2,1]);
    
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

    Mssvec(count_u-MeshU.nDOFe^2:count_u-1) = Mss_e;
    Mffvec(count_u-MeshU.nDOFe^2:count_u-1) = Mff_e;
    Msfvec(count_u-MeshU.nDOFe^2:count_u-1) = Msf_e;

    Cspvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Csp_e;
    Cppvec(count_p-MeshP.nDOFe^2:count_p-1) = Cpp_e;
    
    Kssvec(count_u-MeshU.nDOFe^2:count_u-1) = Kss_e;
    Kspvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Ksp_e;
    Kfpvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Kfp_e;
    Kffvec(count_u-MeshU.nDOFe^2:count_u-1) = Kff_e;

    % u-u
    rowu(count_u-MeshU.nDOFe^2:count_u-1) = rowu_e;
    colu(count_u-MeshU.nDOFe^2:count_u-1) = colu_e;
    % p-p
    rowp(count_p-MeshP.nDOFe^2:count_p-1) = rowp_e;
    colp(count_p-MeshP.nDOFe^2:count_p-1) = colp_e;
    % u-p
    rowup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = rowup_e;
    colup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = colup_e;
end

% sparse matrices
Mss = sparse(rowu, colu, Mssvec, MeshU.nDOF, MeshU.nDOF);
Mff = sparse(rowu, colu, Mffvec, MeshU.nDOF, MeshU.nDOF);
Msf = sparse(rowu, colu, Msfvec, MeshU.nDOF, MeshU.nDOF);
Mfs = Msf.';

Csp = sparse(rowup, colup, Cspvec, MeshU.nDOF, MeshP.nDOF);
Cps = Csp.';
Cpp = sparse(rowp, colp, Cppvec, MeshP.nDOF, MeshP.nDOF);
Cfs = sparse(MeshU.nDOF, MeshU.nDOF);

Kss = sparse(rowu, colu, Kssvec, MeshU.nDOF, MeshU.nDOF);
Ksp = sparse(rowup, colup, Kspvec, MeshU.nDOF, MeshP.nDOF);
Kfp = sparse(rowup, colup, Kfpvec, MeshU.nDOF, MeshP.nDOF);
Kff = sparse(rowu, colu, Kffvec, MeshU.nDOF, MeshU.nDOF);
Kpf = Kfp.';

end