function [Muu, Mpu, Mnu, Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun] = ComputeMatricesDynamic_Spanos(Material, MeshU, MeshP, MeshN, Quad)
% ------------------------------------------------------------------------
% Compute System Matrices for 1D dynamic simulation
% ------------------------------------------------------------------------
% Input parameters: Material, Mesh, Control, Quad
% ------------------------------------------------------------------------
% version 3: correcting matrices for Spanos formulation
% ------------------------------------------------------------------------

ne = MeshU.ne; % number of elements
nq = Quad.nq; % total number of integration points

% constitutive matrix
C = getConstitutiveMatrix(Material, MeshU);

%% Initialize global matrices
% initialize vector sizes
% u-u 
rowu = zeros(ne*MeshU.nDOFe^2,1);
colu = zeros(ne*MeshU.nDOFe^2,1);
% p-p
rowp = zeros(ne*MeshP.nDOFe^2,1);
colp = zeros(ne*MeshP.nDOFe^2,1);
% n-n
rown = zeros(ne*MeshN.nDOFe^2,1);
coln = zeros(ne*MeshN.nDOFe^2,1);
% u-p
rowup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
colup = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
% p-u
rowpu = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
colpu = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
% u-n
rowun = zeros(ne*MeshU.nDOFe*MeshN.nDOFe,1);
colun = zeros(ne*MeshU.nDOFe*MeshN.nDOFe,1);
% n-u
rownu = zeros(ne*MeshU.nDOFe*MeshN.nDOFe,1);
colnu = zeros(ne*MeshU.nDOFe*MeshN.nDOFe,1);
% p-n
rowpn = zeros(ne*MeshP.nDOFe*MeshN.nDOFe,1);
colpn = zeros(ne*MeshP.nDOFe*MeshN.nDOFe,1);
% n-p
rownp = zeros(ne*MeshP.nDOFe*MeshN.nDOFe,1);
colnp = zeros(ne*MeshP.nDOFe*MeshN.nDOFe,1);

% initialize vectorized global matrices
Kuuvec = zeros(ne*MeshU.nDOFe^2,1);
Kppvec = zeros (ne*MeshP.nDOFe^2,1);
Svec = zeros (ne*MeshP.nDOFe^2,1);
Knnvec = zeros (ne*MeshN.nDOFe^2,1);
Kupvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
Kpuvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
Knuvec = zeros(ne*MeshU.nDOFe*MeshN.nDOFe,1);
Kpnvec = zeros(ne*MeshP.nDOFe*MeshN.nDOFe,1);
Knpvec = zeros(ne*MeshP.nDOFe*MeshN.nDOFe,1);
% mass matrices
Muuvec = zeros(ne*MeshU.nDOFe^2,1);
Mpuvec = zeros(ne*MeshU.nDOFe*MeshP.nDOFe,1);
Mnuvec = zeros(ne*MeshU.nDOFe*MeshN.nDOFe,1);

% DOF counter
count_u = 1;
count_p = 1;
count_n = 1;
count_up = 1;
count_un = 1;
count_pn = 1;

%% Coupled matrices
for e = 1:ne
    % element connectivity
    connu_e = MeshU.conn(e,:);
    connu_e = reshape(connu_e',MeshU.nne,[]);
    connp_e = MeshP.conn(e,:);
    connp_e = reshape(connp_e',MeshP.nne,[]);
    connN_e = MeshN.conn(e,:);
    connN_e = reshape(connN_e',MeshN.nne,[]);
    
    % element DOF numbers
    dofu_e = MeshU.DOF(connu_e,:);
    dofu_e = reshape(dofu_e',MeshU.nDOFe,[]);
    dofp_e = MeshP.DOF(connp_e,:);
    dofp_e = reshape(dofp_e',MeshP.nDOFe,[]);
    dofn_e = MeshP.DOF(connN_e,:);
    dofn_e = reshape(dofn_e',MeshN.nDOFe,[]);
    
    % global coordinates
    gcoordsU = MeshU.coords(connu_e,:);
    gcoordsP = MeshP.coords(connp_e,:);
    gcoordsN = MeshN.coords(connN_e,:);

    % initialize local matrices
    Kuu_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Kpp_e = zeros(MeshP.nDOFe, MeshP.nDOFe);
    S_e = zeros(MeshP.nDOFe, MeshP.nDOFe);
    Knn_e = zeros(MeshN.nDOFe, MeshN.nDOFe);
    Kup_e = zeros(MeshU.nDOFe,MeshP.nDOFe);
    Kpu_e = zeros(MeshP.nDOFe, MeshU.nDOFe);
    Kpn_e = zeros(MeshP.nDOFe, MeshN.nDOFe);
    Knp_e = zeros(MeshN.nDOFe, MeshP.nDOFe);
    Knu_e = zeros(MeshN.nDOFe, MeshU.nDOFe);
    Muu_e = zeros(MeshU.nDOFe, MeshU.nDOFe);
    Mpu_e = zeros(MeshP.nDOFe, MeshU.nDOFe);
    Mnu_e = zeros(MeshN.nDOFe, MeshU.nDOFe);
    
    % loop over integration points
    for ip = 1:nq

        % N matrices
        Nu = getN(MeshU, Quad, ip);
        Np = getN(MeshP, Quad, ip);
        Nn = getN(MeshN, Quad, ip);
        
        % N derivatives
        dNu = getdN(MeshU, Quad, ip);
        dNp = getdN(MeshP, Quad, ip);
        dNn = getdN(MeshN, Quad, ip);

        % Jacobian matrix
        Ju = dNu*gcoordsU;
        Jp = dNp*gcoordsP;
        Jn = dNn*gcoordsN;
        % Jacobian determinant
        Jdet = det(Jp);

        % B matrices
        Bu = Ju\dNu;
        Bp = Jp\dNp;
        Bn = Jn\dNn;
        
        % changing matrices to Voigt form
        NuVoigt = getNVoigt(MeshU, Nu);
        NpVoigt = getNVoigt(MeshP, Np);
        NnVoigt = getNVoigt(MeshN, Nn);
        
        BuVoigt = getBVoigt(MeshU, Bu);
        BpVoigt = getBVoigt(MeshP, Bp);      
        BnVoigt = getBVoigt(MeshN, Bn);
        
        % assemble local square matrices
        Kuu_e = Kuu_e + (BuVoigt.') * C * BuVoigt * Material.t * Jdet * Quad.w(ip,1);
        Kpp_e = Kpp_e + Material.kf * (BpVoigt.') * BpVoigt * Material.t * Jdet * Quad.w(ip,1);
        S_e = S_e + (Material.n / Material.Kf) * (NpVoigt.') * NpVoigt * Material.t * Jdet * Quad.w(ip,1);
        Knn_e = Knn_e + (NnVoigt.') * NnVoigt * Material.t * Jdet * Quad.w(ip,1);
        Muu_e = Muu_e + Material.rho * (NuVoigt.') * NuVoigt * Material.t * Jdet * Quad.w(ip,1);
        
        % assemble local square matrices p-n
        Kpn_e = Kpn_e + (NpVoigt.') * NnVoigt  * Material.t * Jdet * Quad.w(ip,1);
        Knp_e = Knp_e + Material.deltaF * Material.kf / Material.n * (BnVoigt.') * BpVoigt * Material.t * Jdet * Quad.w(ip,1);
                
        if MeshU.nsd == 2
            m = [1; 1; 0]; % mapping vector for plane stress
            Kup_e = Kup_e + Material.alpha * (BuVoigt.') * m * NpVoigt * Material.t * Jdet * Quad.w(ip,1);
            Kpu_e = Kpu_e + Material.n * (NpVoigt.') * (m.') * BuVoigt * Material.t * Jdet * Quad.w(ip,1);
            Knu_e = Knu_e + (Material.deltaF - Material.deltaS) * (NnVoigt.') * (m.') * BuVoigt  * Material.t * Jdet * Quad.w(ip,1);
            Mpu_e = Mpu_e + Material.rho_f * Material.kf * (NpVoigt.') * (m.') * BuVoigt * Material.t * Jdet * Quad.w(ip,1);
            Mnu_e = Mnu_e + Material.deltaF * Material.rho_f * Material.kf / Material.n * (NnVoigt.') * (m.') * BuVoigt * Material.t * Jdet * Quad.w(ip,1);
        else
            Kup_e = Kup_e + Material.alpha * (BuVoigt.') * NpVoigt * Material.t * Jdet * Quad.w(ip,1);
            Kpu_e = Kpu_e + Material.n * (NpVoigt.') * BuVoigt * Material.t * Jdet * Quad.w(ip,1);
            Knu_e = Knu_e + (Material.deltaF - Material.deltaS) * (NnVoigt.') * BuVoigt  * Material.t * Jdet * Quad.w(ip,1);
            Mpu_e = Mpu_e + Material.rho_f * Material.kf * (NpVoigt.') * BuVoigt * Material.t * Jdet * Quad.w(ip,1);
            Mnu_e = Mnu_e + Material.deltaF * Material.rho_f * Material.kf / Material.n * (NnVoigt.') * BuVoigt * Material.t * Jdet * Quad.w(ip,1);
        end
    end

    % vectorized matrices
    count_u = count_u + MeshU.nDOFe^2;
    count_p = count_p + MeshP.nDOFe^2;
    count_n = count_n + MeshN.nDOFe^2;
    count_up = count_up + MeshU.nDOFe*MeshP.nDOFe;
    count_un = count_un + MeshU.nDOFe*MeshN.nDOFe;
    count_pn = count_pn + MeshP.nDOFe*MeshN.nDOFe;
    
    Kuu_e = reshape(Kuu_e, [MeshU.nDOFe^2,1]);
    Kpp_e = reshape(Kpp_e, [MeshP.nDOFe^2,1]);
    S_e = reshape(S_e, [MeshP.nDOFe^2,1]);
    Knn_e = reshape(Knn_e, [MeshN.nDOFe^2,1]);
    Kup_e = reshape(Kup_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Kpu_e = reshape(Kpu_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Knu_e = reshape(Knu_e, [MeshU.nDOFe*MeshN.nDOFe,1]);
    Knp_e = reshape(Knp_e, [MeshP.nDOFe*MeshN.nDOFe,1]);
    Kpn_e = reshape(Kpn_e, [MeshP.nDOFe*MeshN.nDOFe,1]);
    Muu_e = reshape(Muu_e, [MeshU.nDOFe^2,1]);
    Mpu_e = reshape(Mpu_e, [MeshU.nDOFe*MeshP.nDOFe,1]);
    Mnu_e = reshape(Mnu_e, [MeshU.nDOFe*MeshN.nDOFe,1]);
    
    % u-u
    rowmatrix_u = dofu_e*ones(1,MeshU.nDOFe);
    rowu_e = reshape(rowmatrix_u, [MeshU.nDOFe^2,1]);
    colu_e = reshape(rowmatrix_u', [MeshU.nDOFe^2,1]);
    % p-p
    rowmatrix_p = dofp_e*ones(1,MeshP.nDOFe);
    rowp_e = reshape(rowmatrix_p, [MeshP.nDOFe^2,1]);
    colp_e = reshape(rowmatrix_p', [MeshP.nDOFe^2,1]);
    % n-n
    rowmatrix_n = dofn_e*ones(1,MeshN.nDOFe);
    rown_e = reshape(rowmatrix_n, [MeshN.nDOFe^2,1]);
    coln_e = reshape(rowmatrix_n', [MeshN.nDOFe^2,1]);
    % u-p
    rowup_e = reshape(dofu_e*ones(1,MeshP.nDOFe),[MeshU.nDOFe*MeshP.nDOFe,1]);
    colup_e = reshape(ones(MeshU.nDOFe,1)*dofp_e.',[MeshU.nDOFe*MeshP.nDOFe,1]);
    % p-u
    rowpu_e = reshape(dofp_e*ones(1,MeshU.nDOFe),[MeshU.nDOFe*MeshP.nDOFe,1]);
    colpu_e = reshape(ones(MeshP.nDOFe,1)*dofu_e.',[MeshU.nDOFe*MeshP.nDOFe,1]);
    % u-n
    rowun_e = reshape(dofu_e*ones(1,MeshN.nDOFe),[MeshU.nDOFe*MeshN.nDOFe,1]);
    colun_e = reshape(ones(MeshU.nDOFe,1)*dofn_e.',[MeshU.nDOFe*MeshN.nDOFe,1]);
    % n-u
    rownu_e = reshape(dofn_e*ones(1,MeshU.nDOFe),[MeshU.nDOFe*MeshN.nDOFe,1]);
    colnu_e = reshape(ones(MeshN.nDOFe,1)*dofu_e.',[MeshU.nDOFe*MeshN.nDOFe,1]);
    % p-n
    rowpn_e = reshape(dofp_e*ones(1,MeshN.nDOFe),[MeshP.nDOFe*MeshN.nDOFe,1]);
    colpn_e = reshape(ones(MeshP.nDOFe,1)*dofn_e.',[MeshP.nDOFe*MeshN.nDOFe,1]);
    % n-p
    rownp_e = reshape(dofn_e*ones(1,MeshP.nDOFe),[MeshP.nDOFe*MeshN.nDOFe,1]);
    colnp_e = reshape(ones(MeshN.nDOFe,1)*dofp_e.',[MeshP.nDOFe*MeshN.nDOFe,1]);

    Kuuvec(count_u-MeshU.nDOFe^2:count_u-1) = Kuu_e;
    Kppvec(count_p-MeshP.nDOFe^2:count_p-1) = Kpp_e;
    Svec(count_p-MeshP.nDOFe^2:count_p-1) = S_e;
    Knnvec(count_n-MeshN.nDOFe^2:count_n-1) = Knn_e;
    Kupvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Kup_e;
    Kpuvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Kpu_e;
    Knuvec(count_un-MeshU.nDOFe*MeshN.nDOFe:count_un-1) = Knu_e;
    Knpvec(count_pn-MeshP.nDOFe*MeshN.nDOFe:count_pn-1) = Knp_e;
    Kpnvec(count_pn-MeshP.nDOFe*MeshN.nDOFe:count_pn-1) = Kpn_e;
    Muuvec(count_u-MeshU.nDOFe^2:count_u-1) = Muu_e;
    Mpuvec(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = Mpu_e;
    Mnuvec(count_un-MeshU.nDOFe*MeshN.nDOFe:count_un-1) = Mnu_e;

    % u-u
    rowu(count_u-MeshU.nDOFe^2:count_u-1) = rowu_e;
    colu(count_u-MeshU.nDOFe^2:count_u-1) = colu_e;
    % p-p
    rowp(count_p-MeshP.nDOFe^2:count_p-1) = rowp_e;
    colp(count_p-MeshP.nDOFe^2:count_p-1) = colp_e;
    % n-n
    rown(count_n-MeshN.nDOFe^2:count_n-1) = rown_e;
    coln(count_n-MeshN.nDOFe^2:count_n-1) = coln_e;
    % u-p
    rowup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = rowup_e;
    colup(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = colup_e;
    % p-u
    rowpu(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = rowpu_e;
    colpu(count_up-MeshU.nDOFe*MeshP.nDOFe:count_up-1) = colpu_e;
    % u-n
    rowun(count_un-MeshU.nDOFe*MeshN.nDOFe:count_un-1) = rowun_e;
    colun(count_un-MeshU.nDOFe*MeshN.nDOFe:count_un-1) = colun_e;
    % n-u
    rownu(count_un-MeshU.nDOFe*MeshN.nDOFe:count_un-1) = rownu_e;
    colnu(count_un-MeshU.nDOFe*MeshN.nDOFe:count_un-1) = colnu_e;
    % p-n
    rowpn(count_pn-MeshP.nDOFe*MeshN.nDOFe:count_pn-1) = rowpn_e;
    colpn(count_pn-MeshP.nDOFe*MeshN.nDOFe:count_pn-1) = colpn_e;
    % n-p
    rownp(count_pn-MeshP.nDOFe*MeshN.nDOFe:count_pn-1) = rownp_e;
    colnp(count_pn-MeshP.nDOFe*MeshN.nDOFe:count_pn-1) = colnp_e;
end

% sparse square matrices
Kuu = sparse(rowu, colu, Kuuvec, MeshU.nDOF, MeshU.nDOF);
Kpp = sparse(rowp, colp, Kppvec, MeshP.nDOF, MeshP.nDOF);
S = sparse(rowp, colp, Svec, MeshP.nDOF, MeshP.nDOF);
Knn = sparse(rown, coln, Knnvec, MeshN.nDOF, MeshN.nDOF);
Kup = sparse(rowup, colup, Kupvec, MeshU.nDOF, MeshP.nDOF);
Kpu = sparse(rowpu, colpu, Kpuvec, MeshP.nDOF, MeshU.nDOF);
Knu = sparse(rownu, colnu, Knuvec, MeshN.nDOF, MeshU.nDOF);
Knp = sparse(rownp, colnp, Knpvec, MeshN.nDOF, MeshP.nDOF);
Kpn = sparse(rowpn, colpn, Kpnvec, MeshP.nDOF, MeshN.nDOF);
Muu = sparse(rowu, colu, Muuvec, MeshU.nDOF, MeshU.nDOF);
Mpu = sparse(rowpu, colpu, Mpuvec, MeshP.nDOF, MeshU.nDOF);
Mnu = sparse(rownu, colnu, Mnuvec, MeshN.nDOF, MeshU.nDOF);
Kun = sparse(MeshU.nDOF, MeshN.nDOF);

end