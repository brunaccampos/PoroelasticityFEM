function ComputeError()
% ------------------------------------------------------------------------
% Compute error using the Le norm of the displacements
% ------------------------------------------------------------------------

clearvars
clear, clc
close all

%% Number of sims
nsims = 5;

%% Mesh 1
% load file
load Results_m1ManufacturedL3_EN.mat
% mesh
MeshP1 = MeshP;
MeshU1 = MeshU;
% FEM approximation
dp1 = p;
du1 = u;
e1 = e;
s1 = s;
q1 = q;
% exact solution
dp_exact1 = p_an;
du_exact1 = u_an;
e_exact1 = e_an;
s_exact1 = s_an;
q_exact1 = q_an;

%% Mesh 2
load Results_m2ManufacturedL3_EN.mat
% mesh
MeshP2 = MeshP;
MeshU2 = MeshU;
% FEM approximation
dp2 = p;
du2 = u;
e2 = e;
s2 = s;
q2 = q;
% exact solution
dp_exact2 = p_an;
du_exact2 = u_an;
e_exact2 = e_an;
s_exact2 = s_an;
q_exact2 = q_an;

%% Mesh 3
load Results_m3ManufacturedL3_EN.mat
% mesh
MeshP3 = MeshP;
MeshU3 = MeshU;
% FEM approximation
dp3 = p;
du3 = u;
e3 = e;
s3 = s;
q3 = q;
% exact solution
dp_exact3 = p_an;
du_exact3 = u_an;
e_exact3 = e_an;
s_exact3 = s_an;
q_exact3 = q_an;

%% Mesh 4
load Results_m4ManufacturedL3_EN.mat
% mesh
MeshP4 = MeshP;
MeshU4 = MeshU;
% FEM approximation
dp4 = p;
du4 = u;
e4 = e;
s4 = s;
q4 = q;
% exact solution
dp_exact4 = p_an;
du_exact4 = u_an;
e_exact4 = e_an;
s_exact4 = s_an;
q_exact4 = q_an;
% 
%% Mesh 5
load Results_m5ManufacturedL3_EN.mat
% mesh
MeshP5 = MeshP;
MeshU5 = MeshU;
% FEM approximation
dp5 = p;
du5 = u;
e5 = e;
s5 = s;
q5 = q;
% exact solution
dp_exact5 = p_an;
du_exact5 = u_an;
e_exact5 = e_an;
s_exact5 = s_an;
q_exact5 = q_an;

% [p_an, u_an] = getAnalyticResult_Symb(Material, MeshU, MeshP, BC, Control);


%% Compute L2 error norm
% initialize variables
eL2p = zeros(nsims,1);
eL2u = zeros(nsims,1);
eENp = zeros(nsims,1);
eENu = zeros(nsims,1);
h = zeros(nsims,1);
normu = zeros(nsims,1);
normp = zeros(nsims,1);

% loop over meshes
for sim = 1:nsims
    
    switch sim
        case 1
            dp = dp1;
            du = du1;
            e = e1;
            s = s1;
            q = q1;
            dp_exact = dp_exact1;
            du_exact = du_exact1;
            e_exact = e_exact1;
            s_exact = s_exact1;
            q_exact = q_exact1;
            MeshP = MeshP1;
            MeshU = MeshU1;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,'inf');
        case 2
            dp = dp2;
            du = du2;
            e = e2;
            s = s2;
            q = q2;
            dp_exact = dp_exact2;
            du_exact = du_exact2;
            e_exact = e_exact2;
            s_exact = s_exact2;
            q_exact = q_exact2;
            MeshP = MeshP2;
            MeshU = MeshU2;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,'inf');
        case 3
            dp = dp3;
            du = du3;
            e = e3;
            s = s3;
            q = q3;
            dp_exact = dp_exact3;
            du_exact = du_exact3;
            e_exact = e_exact3;
            s_exact = s_exact3;
            q_exact = q_exact3;
            MeshP = MeshP3;
            MeshU = MeshU3;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,'inf');
        case 4
            dp = dp4;
            du = du4;
            e = e4;
            s = s4;
            q = q4;
            dp_exact = dp_exact4;
            du_exact = du_exact4;
            e_exact = e_exact4;
            s_exact = s_exact4;
            q_exact = q_exact4;
            MeshP = MeshP4;
            MeshU = MeshU4;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,'inf');
        case 5
            dp = dp5;
            du = du5;
            e = e5;
            s = s5;
            q = q5;
            dp_exact = dp_exact5;
            du_exact = du_exact5;
            e_exact = e_exact5;
            s_exact = s_exact5;
            q_exact = q_exact5;
            MeshP = MeshP5;
            MeshU = MeshU5;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,'inf');
        case 6
            dp = dp6;
            du = du6;
            e = e6;
            s = s6;
            q = q6;
            dp_exact = dp_exact6;
            du_exact = dp_exact6;
            e_exact = e_exact6;
            s_exact = s_exact6;
            q_exact = q_exact6;
            MeshP = MeshP6;
            MeshU = MeshU6;
            normu(sim,1) = norm(du-du_exact,2);
            normp(sim,1) = norm(dp-dp_exact,'inf');
    end
    
    % Mesh size
    switch MeshP.nsd
        case 1
            h(sim) = max(MeshP.coords)/MeshP.ne;
            MeshP.xdofs = MeshP.DOF;
            MeshP.ydofs = zeros(length(MeshP.DOF),1);
            MeshU.xdofs = MeshU.DOF;
            MeshU.ydofs = zeros(length(MeshU.DOF),1);
        case 2
            gcoordsp = MeshP.coords(MeshP.conn(1,:),:);
            h(sim) = sqrt(polyarea(gcoordsp(:,1),gcoordsp(:,2)));
    end
    
%     % exact solution
    d_exact = @(x) x.^5 - x.^4;
    e_exact = @(x) 5*x.^4 - 4*x.^3;
    s_exact = e_exact;
    
    Control.nqU = 8;
    MeshU.field = 'u';
    Quad = GlobalQuad(MeshU, Control);
    
    % Calculate error norms
    % L2 norm
    eL2p_num = 0;
    eL2p_den = 0;
    eL2u_num = 0;
    eL2u_den = 0;
    % energy norm
    eENp_num = 0;
    eENp_den = 0;
    eENu_num = 0;
    eENu_den = 0;
    
    % loop over elements
    for i = 1:MeshP.ne
        % element connectivity
        connp_e = MeshP.conn(i,:);
        connu_e = MeshU.conn(i,:);
        % global coordinates
        gcoordsp = MeshP.coords(connp_e,:);
                
        % loop over integration points
        for ip = 1:Quad.nq
            % N matrices
            Np = getN(MeshP, Quad, ip);
            Nu = getN(MeshU, Quad, ip);
            % N derivatives
            dNp = getdN(MeshP, Quad, ip);
            % Jacobian matrix
            J = dNp*gcoordsp;
            % Jacobian determinant
            Jdet = det(J);
            
            % approximated displacement at quadrature point
            switch MeshP.nsd
                case 1
                    % approximated displacement at quadrature point
                    uxh_p = Np*dp(MeshP.xdofs(connp_e)');
                    uyh_p = zeros(length(uxh_p));
                    uxh_u = Nu*du(MeshU.xdofs(connu_e)');
                    uyh_u = zeros(length(uxh_u));
                    % exact displacement at quadrature point
                    uxe_p = Np*dp_exact(MeshP.xdofs(connp_e));
                    uye_p = zeros(length(uxe_p));
                    uxe_u = Nu*du_exact(MeshU.xdofs(connu_e));
                    uye_u = zeros(length(uxe_u));

%                     uxe_p = eval(subs(p_an, Np*MeshP.coords(connp_e)));
%                     uxe_u = eval(subs(u_an, Nu*MeshU.coords(connu_e)));

                    uxe_u = eval(subs(d_exact, Nu*MeshU.coords(connu_e)));

                case 2
                    % approximated displacement at quadrature point
                    uxh_p = Np*dp(MeshP.xdofs(connp_e)');
                    uyh_p = Np*dp(MeshP.ydofs(connp_e)');
                    uxh_u = Nu*du(MeshU.xdofs(connu_e)');
                    uyh_u = Nu*du(MeshU.ydofs(connu_e)');
                    % exact displacement at quadrature point
                    uxe_p = Np*dp_exact(MeshP.xdofs(connp_e)');
                    uye_p = Np*dp_exact(MeshP.ydofs(connp_e)');
                    uxe_u = Nu*du_exact(MeshU.xdofs(connu_e)');
                    uye_u = Nu*du_exact(MeshU.ydofs(connu_e)');
            end

            % L2 norm pressure
            eL2p_num = eL2p_num + [uxh_p - uxe_p , uyh_p - uye_p] * [uxh_p - uxe_p ; uyh_p - uye_p] * Quad.w(ip,1) * Jdet;
            eL2p_den = eL2p_den + [uxe_p , uye_p] * [uxe_p ; uye_p] * Quad.w(ip,1) * Jdet;
            % L2 norm displacement
            eL2u_num = eL2u_num + [uxh_u - uxe_u , uyh_u - uye_u] * [uxh_u - uxe_u ; uyh_u - uye_u] * Quad.w(ip,1) * Jdet;
            eL2u_den = eL2u_den + [uxe_u , uye_u] * [uxe_u ; uye_u] * Quad.w(ip,1) * Jdet;
            
            % approximated strain, stress, and flux at quadrature points
            eh = Nu*e(connu_e);
            sh = Nu*s(connu_e);
            qh = Np*q(connp_e);
            
            % exact strain, stress, and flux at quadrature points
%             ee = Nu*e_exact(connu_e);
%             se = Nu*s_exact(connu_e);
            qe = Np*q_exact(connp_e);
            
            ee = eval(subs(e_exact, Nu*MeshU.coords(connu_e)));
            se = eval(subs(s_exact, Nu*MeshU.coords(connu_e)));

            % energy norm
            eENu_num = eENu_num + (eh-ee) * (sh-se) * Quad.w(ip,1) * Jdet;
            eENu_den = eENu_den + ee * se * Quad.w(ip,1) * Jdet;
            eENp_num = eENp_num + (qh-qe) * (qh-qe) * Quad.w(ip,1) * Jdet;
            eENp_den = eENp_den + qe * qe * Quad.w(ip,1) * Jdet;
        end
    end
    
    eL2p(sim) = sqrt(eL2p_num/eL2p_den);
    eL2u(sim) = sqrt(eL2u_num/eL2u_den);    
    
    eENp(sim) = sqrt(eENp_num/eENp_den);    
    eENu(sim) = sqrt(eENu_num/eENu_den);
end

% Determine slope of L2 norm
pL2p = polyfit(log(h), log(eL2p),1);
m_L2p = pL2p(1);

pL2u = polyfit(log(h), log(eL2u),1);
m_L2u = pL2u(1);

pENp = polyfit(log(h), log(eENp),1);
m_ENp = pENp(1);

pENu = polyfit(log(h), log(eENu),1);
m_ENu = pENu(1);


%% Step 2 - Calculate the slope of each curve
% L2 norm pressure
figure;
loglog(h,eL2p,'m-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('%s - %s Convergence for pressure order %.2f', MeshP.type, MeshP.field, m_L2p));
% L2 norm displacement
figure;
loglog(h,eL2u,'g-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('L2-norm');
title(sprintf('%s - %s Convergence for displacement order %.2f', MeshU.type, MeshU.field, m_L2u));

% Energy norm pressure
figure;
loglog(h,eENp,'b-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('Energy norm');
title(sprintf('%s - %s Convergence of EN-norm: %.2f', MeshP.type, MeshP.field, m_ENp));
% Energy norm displacement
figure;
loglog(h,eENu,'k-o', 'LineWidth', 1.5);
xlabel('Mesh size (m)');
ylabel('Energy norm');
title(sprintf('%s - %s Convergence of EN-norm: %.2f', MeshU.type, MeshU.field, m_ENu));


end