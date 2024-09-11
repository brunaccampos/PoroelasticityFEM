function [ErrorComp] = ComputeMeshSizeError_UPW(MeshU, MeshP, Solution, Plot, Control)
% ------------------------------------------------------------------------
% Compute solid displacement, fluid displacement and pressure errors 
% related to the mesh size in 1D
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
ds = Solution.u;
dp = Solution.p;
df = Solution.w;

% exact solutions evaluated at nodes
ds_exact = Plot.uan_space;
dp_exact = Plot.pan_space;
df_exact = Plot.wan_space;

% exact solution with symbolic function
s_exact = @(x) Control.uan_symb(x,Control.t);
p_exact = @(x) Control.pan_symb(x,Control.t);
f_exact = @(x) Control.wan_symb(x,Control.t);

% higher quadrature to evaluate error
Control.nqU = 8;
Quad = GlobalQuad(MeshU, Control);

% numerator and denominator terms
eL2us_num = 0;
eL2us_den = 0;
eL2p_num = 0;
eL2p_den = 0;
eL2w_num = 0;
eL2w_den = 0;

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
       ush = Nu*ds(MeshU.DOF(connu_e)');
       ph = Np*dp(MeshP.DOF(connp_e)');
       ufh = Nu*df(MeshU.DOF(connu_e)');

       % exact
       use = Nu*ds_exact(MeshU.DOF(connu_e)');
       pe = Np*dp_exact(MeshP.DOF(connp_e)');
       ufe = Nu*df_exact(MeshU.DOF(connu_e)');

       % exact with function handle
%        use = eval(subs(s_exact, Nu*MeshU.coords(connu_e)));
%        pe = eval(subs(p_exact, Np*MeshP.coords(connp_e)));
%        ufe = eval(subs(f_exact, Nu*MeshU.coords(connu_e)));

       % L2 norm solid displacement
       eL2us_num = eL2us_num + (ush - use)^2 * Quad.w(ip,1) * Jdet;
       eL2us_den = eL2us_den + (use)^2 * Quad.w(ip,1) * Jdet;
       
       % L2 norm pressure
       eL2p_num = eL2p_num + (ph - pe)^2 * Quad.w(ip,1) * Jdet;
       eL2p_den = eL2p_den + (pe)^2 * Quad.w(ip,1) * Jdet;
       
       % L2 norm solid displacement
       eL2w_num = eL2w_num + (ufh - ufe)^2 * Quad.w(ip,1) * Jdet;
       eL2w_den = eL2w_den + (ufe)^2 * Quad.w(ip,1) * Jdet;
    end
end

% square root
ErrorComp.eL2u = sqrt(eL2us_num/eL2us_den);
ErrorComp.eL2p = sqrt(eL2p_num/eL2p_den);
ErrorComp.eL2w = sqrt(eL2w_num/eL2w_den);

end