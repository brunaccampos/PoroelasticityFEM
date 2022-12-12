function [Kuu, Kup, Kpp, Kpu, S, Kpn, Knn, Knu, Knp, Kun] = ComputeSystemMatrices_Steady_v2(Material, MeshU, MeshP, MeshN, Control, Quad)
% Compute System Matrices for 1D quasi-steady simulation
% Input parameters: Material, Mesh, Control, Quad
% Output matrices: Kuu, Kup, Kpp, S

ne = MeshU.ne; % number of elements
nq = Control.nq^MeshU.nsd; % total number of integration points

% constitutive matrix
C = getConstitutiveMatrix(Material, MeshU);

%% Initialize global matrices
% initialize vector sizes
rowu = zeros(ne*MeshU.nDOFe^2,1);
colu = zeros(ne*MeshU.nDOFe^2,1);

rowp = zeros(ne*MeshP.nDOFe^2,1);
colp = zeros(ne*MeshP.nDOFe^2,1);

rown = zeros(ne*MeshN.nDOFe^2,1);
coln = zeros(ne*MeshN.nDOFe^2,1);

Kup = zeros(MeshU.nDOF,MeshP.nDOF);
Kuuvec = zeros(ne*MeshU.nDOFe^2,1);
Kppvec = zeros (ne*MeshP.nDOFe^2,1);
Svec = zeros (ne*MeshP.nDOFe^2,1);

%%% new
Knnvec = zeros (ne*MeshN.nDOFe^2,1);
Kpu = zeros(MeshP.nDOF,MeshU.nDOF);
Kpn = zeros(MeshP.nDOF,MeshN.nDOF);
Knu = zeros(MeshN.nDOF,MeshU.nDOF);
Knp = zeros(MeshN.nDOF,MeshP.nDOF);
Kun = zeros(MeshU.nDOF,MeshN.nDOF);
%%% end

% DOF counter
count_u = 1;
count_p = 1;
count_n = 1;

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
    Kup_e = zeros(MeshU.nDOFe,MeshP.nDOFe);
    Kpp_e = zeros (MeshP.nDOFe, MeshP.nDOFe);
    S_e = zeros (MeshP.nDOFe, MeshP.nDOFe);
    
    %%% new
    Knn_e = zeros (MeshN.nDOFe, MeshN.nDOFe);
    Kpn_e = zeros (MeshP.nDOFe, MeshN.nDOFe);
    Knu_e = zeros (MeshN.nDOFe, MeshU.nDOFe);
    Knp_e = zeros (MeshN.nDOFe, MeshP.nDOFe);
    Kpu_e = zeros(MeshP.nDOFe,MeshU.nDOFe);
    %%% end
    
    % loop over integration points
    for ip = 1:nq

        % N matrices
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
        NpVoigt = getNVoigt(MeshP, Np);
        BuVoigt = getBVoigt(MeshU, Bu);
        BpVoigt = getBVoigt(MeshP, Bp);

        NnVoigt = getNVoigt(MeshN, Nn);
        BnVoigt = getBVoigt(MeshN, Bn);
        
        % assemble local matrices
        Kuu_e = Kuu_e + (BuVoigt.') * C * BuVoigt * Jdet * Quad.w(ip,1);
        Kpp_e = Kpp_e + Material.kf * (BpVoigt.') * BpVoigt * Jdet * Quad.w(ip,1);
%         S_e = S_e + Material.Qinv * (NpVoigt.') * NpVoigt * Jdet * Quad.w(ip,1);
        S_e = S_e + (Material.n / Material.Kf) * (NpVoigt.') * NpVoigt * Jdet * Quad.w(ip,1);
        
        %%% new
        Knn_e = Knn_e + (NnVoigt.') * NnVoigt * Jdet * Quad.w(ip,1);
        %%%% WORKING ONLY FOR 1D FOR NOW
        Kpn_e = Kpn_e + (NpVoigt.') * NnVoigt  * Jdet * Quad.w(ip,1);
        Knu_e = Knu_e + (Material.deltaF - Material.deltaS) * (NnVoigt.') * BuVoigt  * Jdet * Quad.w(ip,1);
        Knp_e = Knp_e + Material.deltaF * Material.kf / Material.n * (BnVoigt.') * BpVoigt * Jdet * Quad.w(ip,1);
        %%% end
        
        if MeshU.nsd == 2
            m = [1; 1; 0]; % mapping vector for plane stress
            Kup_e = Kup_e + Material.alpha * (BuVoigt.') * m * NpVoigt * Jdet * Quad.w(ip,1);
        else
            Kup_e = Kup_e + Material.alpha * (BuVoigt.') * NpVoigt * Jdet * Quad.w(ip,1);
            %%% new
            Kpu_e = Kpu_e + Material.n * (NpVoigt.') * BuVoigt * Jdet * Quad.w(ip,1);
            %%% end
        end
    end

    % vectorized matrices
    count_u = count_u + MeshU.nDOFe^2;
    count_p = count_p + MeshP.nDOFe^2;
    count_n = count_n + MeshN.nDOFe^2;
    
    Kuu_e = reshape(Kuu_e, [MeshU.nDOFe^2,1]);
    Kpp_e = reshape(Kpp_e, [MeshP.nDOFe^2,1]);
    S_e = reshape(S_e, [MeshP.nDOFe^2,1]);
    Knn_e = reshape(Knn_e, [MeshN.nDOFe^2,1]);
    
    rowmatrix_u = dofu_e*ones(1,MeshU.nDOFe);
    rowu_e = reshape(rowmatrix_u, [MeshU.nDOFe^2,1]);
    colu_e = reshape(rowmatrix_u', [MeshU.nDOFe^2,1]);

    rowmatrix_p = dofp_e*ones(1,MeshP.nDOFe);
    rowp_e = reshape(rowmatrix_p, [MeshP.nDOFe^2,1]);
    colp_e = reshape(rowmatrix_p', [MeshP.nDOFe^2,1]);

    rowmatrix_n = dofn_e*ones(1,MeshN.nDOFe);
    rown_e = reshape(rowmatrix_n, [MeshN.nDOFe^2,1]);
    coln_e = reshape(rowmatrix_n', [MeshN.nDOFe^2,1]);

    Kuuvec(count_u-MeshU.nDOFe^2:count_u-1) = Kuu_e;
    Kppvec(count_p-MeshP.nDOFe^2:count_p-1) = Kpp_e;
    Svec(count_p-MeshP.nDOFe^2:count_p-1) = S_e;
    Knnvec(count_n-MeshN.nDOFe^2:count_n-1) = Knn_e;

    rowu(count_u-MeshU.nDOFe^2:count_u-1) = rowu_e;
    colu(count_u-MeshU.nDOFe^2:count_u-1) = colu_e;

    rowp(count_p-MeshP.nDOFe^2:count_p-1) = rowp_e;
    colp(count_p-MeshP.nDOFe^2:count_p-1) = colp_e;

    rown(count_n-MeshN.nDOFe^2:count_n-1) = rown_e;
    coln(count_n-MeshN.nDOFe^2:count_n-1) = coln_e;

    % gather matrices
    Lu = getGatherMatrix(MeshU.nDOFe, MeshU.nDOF, e);
    Lp = getGatherMatrix(MeshP.nDOFe, MeshP.nDOF, e);
    Ln = getGatherMatrix(MeshN.nDOFe, MeshN.nDOF, e);
    
    % assemble Kup
    Kup = Kup + (Lu.') * Kup_e * Lp;
    Kpu = Kpu + (Lp.') * Kpu_e * Lu;
    Kpn = Kpn + (Lp.') * Kpn_e * Ln;
    Knu = Knu + (Ln.') * Knu_e * Lu;
    Knp = Knp + (Ln.') * Knp_e * Lp;
end

% sparse matrices
Kuu = sparse(rowu, colu, Kuuvec, MeshU.nDOF, MeshU.nDOF);
Kup = sparse(Kup);
Kpp = sparse(rowp, colp, Kppvec, MeshP.nDOF, MeshP.nDOF);
S = sparse(rowp, colp, Svec, MeshP.nDOF, MeshP.nDOF);

Knn = sparse(rown, coln, Knnvec, MeshN.nDOF, MeshN.nDOF);
Kpn = sparse(Kpn);
Knu = sparse(Knu);
Knp = sparse(Knp);
Kpu = sparse(Kpu);

end