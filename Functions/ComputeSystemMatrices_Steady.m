function [Kuu, Kup, Kpp, S] = ComputeSystemMatrices_Steady(Material, Mesh, Control, Quad)
% Compute System Matrices for 1D quasi-steady simulation
% Input parameters: Material, Mesh, Control, Quad
% Output matrices: Kuu, Kup, Kpp, S

%% Initialize global matrices

% Jacobian
J = Mesh.h/2;

% vector sizes
rowu = zeros(Mesh.ne*Mesh.ndof_u_e^2,1);
colu = zeros(Mesh.ne*Mesh.ndof_u_e^2,1);

rowp = zeros(Mesh.ne*Mesh.ndof_p_e^2,1);
colp = zeros(Mesh.ne*Mesh.ndof_p_e^2,1);

Kup = zeros(Mesh.ndof_u,Mesh.ndof_p);
Kuuvec = zeros(Mesh.ne*Mesh.ndof_u_e^2,1);
Kppvec = zeros (Mesh.ne*Mesh.ndof_p_e^2,1);
Svec = zeros (Mesh.ne*Mesh.ndof_p_e^2,1);

% DOF counter
count_u = 1;
count_p = 1;

%% Coupled matrices
for e = 1:Mesh.ne
    % element connectivity
    connu_e = Mesh.connu(e,:);
    connu_e = reshape(connu_e',Mesh.ndof_u_e,[]);
    connp_e = Mesh.connp(e,:);
    connp_e = reshape(connp_e',Mesh.ndof_p_e,[]);

    % initialize local matrices
    Kuu_e = zeros(Mesh.ndof_u_e, Mesh.ndof_u_e);
    Kup_e = zeros(Mesh.ndof_u_e,Mesh.ndof_p_e);
    Kpp_e = zeros (Mesh.ndof_p_e, Mesh.ndof_p_e);
    S_e = zeros (Mesh.ndof_p_e, Mesh.ndof_p_e);

    % loop over integration points
    for ip = 1:Control.nq

        % N matrices
        Np = getN(Mesh.ndof_p_e, Quad, ip);

        % B matrices
        Bu = getB(Mesh.ndof_u_e, Mesh, Quad, ip);
        Bp = getB(Mesh.ndof_p_e, Mesh, Quad, ip);

        % assemble local matrices
        Kuu_e = Kuu_e + (Bu.') * Material.E * Bu * J * Quad.w_csi(ip,1);
        Kup_e = Kup_e + Material.alpha * (Bu.') * Np * J * Quad.w_csi(ip,1);
        Kpp_e = Kpp_e + Material.kf * (Bp.') * Bp * J * Quad.w_csi(ip,1);
        S_e = S_e + Material.Q * (Np.') * Np * J * Quad.w_csi(ip,1);
    end

    % vectorized matrices
    count_u = count_u + Mesh.ndof_u_e^2;
    count_p = count_p + Mesh.ndof_p_e^2;

    Kuu_e = reshape(Kuu_e, [Mesh.ndof_u_e^2,1]);
    Kpp_e = reshape(Kpp_e, [Mesh.ndof_p_e^2,1]);
    S_e = reshape(S_e, [Mesh.ndof_p_e^2,1]);

    rowmatrix_u = connu_e*ones(1,Mesh.ndof_u_e);
    rowu_e = reshape(rowmatrix_u, [Mesh.ndof_u_e^2,1]);
    colu_e = reshape(rowmatrix_u', [Mesh.ndof_u_e^2,1]);

    rowmatrix_p = connp_e*ones(1,Mesh.ndof_p_e);
    rowp_e = reshape(rowmatrix_p, [Mesh.ndof_p_e^2,1]);
    colp_e = reshape(rowmatrix_p', [Mesh.ndof_p_e^2,1]);

    Kuuvec(count_u-Mesh.ndof_u_e^2:count_u-1) = Kuu_e;
    Kppvec(count_p-Mesh.ndof_p_e^2:count_p-1) = Kpp_e;
    Svec(count_p-Mesh.ndof_p_e^2:count_p-1) = S_e;

    rowu(count_u-Mesh.ndof_u_e^2:count_u-1) = rowu_e;
    colu(count_u-Mesh.ndof_u_e^2:count_u-1) = colu_e;

    rowp(count_p-Mesh.ndof_p_e^2:count_p-1) = rowp_e;
    colp(count_p-Mesh.ndof_p_e^2:count_p-1) = colp_e;

    % gather matrices
    Lu = getGatherMatrix(Mesh.ndof_u_e, Mesh.ndof_u, e);
    Lp = getGatherMatrix(Mesh.ndof_p_e, Mesh.ndof_p, e);

    % assemble Kup
    Kup = Kup + (Lu.') * Kup_e * Lp;
end

% sparse matrices
Kuu = sparse(rowu, colu, Kuuvec, Mesh.ndof_u, Mesh.ndof_u);
Kup = sparse(Kup);
Kpp = sparse(rowp, colp, Kppvec, Mesh.ndof_p, Mesh.ndof_p);
S = sparse(rowp, colp, Svec, Mesh.ndof_p, Mesh.ndof_p);

end