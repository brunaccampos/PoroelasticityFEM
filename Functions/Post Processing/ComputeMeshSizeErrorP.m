function [ErrorP] = ComputeMeshSizeErrorP(MeshP, Solution, Plot, Control)
% ------------------------------------------------------------------------
% Compute error norms
% L2-norm: uses pressure
% Energy norm: uses flux
% H1-norm: uses pressure and flux
% ------------------------------------------------------------------------

%% Initialize variables
% number of elements
ne = MeshP.ne;

% mesh element size
ErrorP.h = max(MeshP.coords)/MeshP.ne;

% approximate solutions
dp = Solution.p;
q = Solution.q;

% exact solutions evaluated at nodes
dp_exact = Plot.pan_space;
q_exact = Solution.q_an;

% exact solution with symbolic function
d_exact = Control.pan_symb;

% higher quadrature to evaluate error
Control.nqP = 16;
Quad = GlobalQuad(MeshP, Control);

% numerator and denominator terms
eL2p_num = 0;
eL2p_den = 0;

eENp_num = 0;
eENp_den = 0;

eH1p_num = 0;
eH1p_den = 0;

%% Compute error norm
% loop over elements
for i = 1:ne
    % element connectivity
    connp_e = MeshP.conn(i,:);
    % global coordinates
    gcoords_p = MeshP.coords(connp_e,:);
    
    % loop over integration points
    for ip = 1:Quad.nq
       % N matrices
       Np = getN(MeshP, Quad, ip);
       % N derivatives
       dNp = getdN(MeshP, Quad, ip);
       % Jacobian matrix
       J = dNp*gcoords_p;
       % Jacobian determinant
       Jdet = det(J);
       
       % values at the quadrature point
       % approximate (FEM)
       pxh = Np*dp(MeshP.DOF(connp_e)');
       % exact
       pxe = Np*dp_exact(MeshP.DOF(connp_e)');
       
       % with function handle
       pxe = eval(subs(d_exact, Np*MeshP.coords(connp_e)));
       
       % L2 norm pressure
       eL2p_num = eL2p_num + (pxh - pxe)^2 * Quad.w(ip,1) * Jdet;
       eL2p_den = eL2p_den + (pxe)^2 * Quad.w(ip,1) * Jdet;
       
       % approximate (FEM)
       qh = Np*q(connp_e);
       
       % exact
       qe = Np*q_exact(connp_e);

%        qe = eval(subs(qex_exact, Np*MeshP.coords(connp_e)));
     
       % energy norm
       eENp_num = eENp_num + (qh-qe) * (qh-qe) * Quad.w(ip,1) * Jdet;
       eENp_den = eENp_den + qe * qe * Quad.w(ip,1) * Jdet;

       % H1 norm
       eH1p_num = eH1p_num + ((qh-qe)^2 + (pxh-pxe)^2) * Quad.w(ip,1) * Jdet;
       eH1p_den = eH1p_den + ((qe)^2 + (pxe)^2) * Quad.w(ip,1) * Jdet;
    end
end

% square root
ErrorP.eL2p = sqrt(eL2p_num/eL2p_den);
ErrorP.eENp = sqrt(eENp_num/eENp_den);    
ErrorP.eH1p = sqrt(eH1p_num/eH1p_den);    

end