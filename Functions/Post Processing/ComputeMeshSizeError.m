function [ErrorComp] = ComputeMeshSizeError(MeshU, MeshP, Solution, Plot, Control)
% ------------------------------------------------------------------------
% Compute displacement and pressure errors related to the mesh size in 1D
% L2-norm: uses displacement / pressure
% Energy norm: uses strain and stress / flux
% H1-norm: uses displacement and strain / pressure and flux
% ------------------------------------------------------------------------

%% Initialize variables
% number of elements
ne = MeshU.ne;

% mesh element size
ErrorComp.h = max(MeshP.coords)/MeshP.ne;

% approximate solutions
du = Solution.u;
dp = Solution.p;
e = Solution.e;
q = Solution.q;

% exact solutions evaluated at nodes
du_exact = Plot.uan_space;
dp_exact = Plot.pan_space;
e_exact = Solution.e_an;
q_exact = Solution.q_an;

% exact solution with symbolic function
d_exact = @(x) Control.uan_symb(x,Control.t);
p_exact = @(x) Control.pan_symb(x,Control.t);

% higher quadrature to evaluate error
Control.nqU = 8;
Quad = GlobalQuad(MeshU, Control);

% numerator and denominator terms
eL2u_num = 0;
eL2u_den = 0;
eL2p_num = 0;
eL2p_den = 0;

eENu_num = 0;
eENu_den = 0;
eENp_num = 0;
eENp_den = 0;

eH1u_num = 0;
eH1u_den = 0;
eH1p_num = 0;
eH1p_den = 0;

%% Compute error norm
% loop over elements
for i = 1:ne
    % element connectivity
    connu_e = MeshU.conn(i,:);
    connp_e = MeshP.conn(i,:);
    % global coordinates
    gcoords_p = MeshP.coords(connp_e,:);
    
    % loop over integration points
    for ip = 1:Quad.nq
       % N matrices
       Nu = getN(MeshU, Quad, ip);
       Np = getN(MeshP, Quad, ip);
       % N derivatives
       dNp = getdN(MeshP, Quad, ip);
       % Jacobian matrix
       J = dNp*gcoords_p;
       % Jacobian determinant
       Jdet = det(J);
       
       % values at the quadrature point
       % approximate (FEM)
       uxh = Nu*du(MeshU.DOF(connu_e)');
       pxh = Np*dp(MeshP.DOF(connp_e)');
       
       % exact
       uxe = Nu*du_exact(MeshU.DOF(connu_e)');
       pxe = Np*dp_exact(MeshP.DOF(connp_e)');
       
       % exact with function handle
       uxe = eval(subs(d_exact, Nu*MeshU.coords(connu_e)));
       pxe = eval(subs(p_exact, Np*MeshP.coords(connp_e)));
       
       % L2 norm displacement
       eL2u_num = eL2u_num + (uxh - uxe)^2 * Quad.w(ip,1) * Jdet;
       eL2u_den = eL2u_den + (uxe)^2 * Quad.w(ip,1) * Jdet;
       
       % L2 norm pressure
       eL2p_num = eL2p_num + (pxh - pxe)^2 * Quad.w(ip,1) * Jdet;
       eL2p_den = eL2p_den + (pxe)^2 * Quad.w(ip,1) * Jdet;
       
       % approximate (FEM)
       eh = Nu*e(connu_e);
       qh = Np*q(connp_e);
       
       % exact
       ee = Nu*e_exact(connu_e);
       qe = Np*q_exact(connp_e);
       
%        ee = eval(subs(e_exact, Nu*MeshU.coords(connu_e)));
%        qe = eval(subs(qex_exact, Np*MeshP.coords(connp_e)));
     
       % energy norm
       eENu_num = eENu_num + (eh-ee)^2 * Quad.w(ip,1) * Jdet;
       eENu_den = eENu_den + (ee)^2 * Quad.w(ip,1) * Jdet;
       eENp_num = eENp_num + (qh-qe)^2 * Quad.w(ip,1) * Jdet;
       eENp_den = eENp_den + (qe)^2 * Quad.w(ip,1) * Jdet;


       % H1 norm
       eH1u_num = eH1u_num + ((eh-ee)^2 + (uxh-uxe)^2) * Quad.w(ip,1) * Jdet;
       eH1u_den = eH1u_den + ((ee)^2 + (uxe)^2) * Quad.w(ip,1) * Jdet;
       eH1p_num = eH1p_num + ((qh-qe)^2 + (pxh-pxe)^2) * Quad.w(ip,1) * Jdet;
       eH1p_den = eH1p_den + ((qe)^2 + (pxe)^2) * Quad.w(ip,1) * Jdet;
    end
end

% square root
ErrorComp.eL2u = sqrt(eL2u_num/eL2u_den);
ErrorComp.eL2p = sqrt(eL2p_num/eL2p_den);

ErrorComp.eENu = sqrt(eENu_num/eENu_den);
ErrorComp.eENp = sqrt(eENp_num/eENp_den);    

ErrorComp.eH1u = sqrt(eH1u_num/eH1u_den);
ErrorComp.eH1p = sqrt(eH1p_num/eH1p_den);    

end