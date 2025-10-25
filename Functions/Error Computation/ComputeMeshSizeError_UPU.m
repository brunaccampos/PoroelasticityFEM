% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [ErrorComp] = ComputeMeshSizeError_UPU(MeshU, MeshP, Solution, Plot, Control)
% Compute solid displacement, fluid displacement and pressure errors 
% related to the mesh size in 1D
% L2-norm: uses displacement / pressure
% Energy norm: uses strain and stress / flux
% H1-norm: uses displacement and strain / pressure and flux

%% Initialize variables
% number of elements
ne = MeshU.ne;

% mesh element size
ErrorComp.h = max(MeshP.coords)/MeshP.ne;

% approximate solutions
ds = Solution.u;
dp = Solution.p;
df = Solution.uf;

strains = Solution.strain;
gradp = Solution.gradp;
strainf = Solution.strainf;

% exact solutions evaluated at nodes
ds_exact = Plot.uan_space;
dp_exact = Plot.pan_space;
df_exact = Plot.ufan_space;

strains_exact = Solution.strain_an;
gradp_exact = Solution.gradp_an;
strainf_exact = Solution.strainf_an;

% exact solution with symbolic function
s_exact = @(x) Control.uan_symb(x,Control.t);
p_exact = @(x) Control.pan_symb(x,Control.t);
f_exact = @(x) Control.ufan_symb(x,Control.t);

% higher quadrature to evaluate error
Control.nqU = 8;
Quad = GlobalQuad(MeshU, Control);

% numerator and denominator terms
eL2us_num = 0;
eL2us_den = 0;
eL2p_num = 0;
eL2p_den = 0;
eL2uf_num = 0;
eL2uf_den = 0;

eENus_num = 0;
eENus_den = 0;
eENp_num = 0;
eENp_den = 0;
eENuf_num = 0;
eENuf_den = 0;

eH1us_num = 0;
eH1us_den = 0;
eH1p_num = 0;
eH1p_den = 0;
eH1uf_num = 0;
eH1uf_den = 0;

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
       eL2uf_num = eL2uf_num + (ufh - ufe)^2 * Quad.w(ip,1) * Jdet;
       eL2uf_den = eL2uf_den + (ufe)^2 * Quad.w(ip,1) * Jdet;

       % approximate (FEM)
       esh = Nu*strains(connu_e);
       gph = Np*gradp(connp_e);
       efh = Nu*strainf(connu_e);

       % exact
       ese = Nu*strains_exact(connu_e);
       gpe = Np*gradp_exact(connp_e);
       efe = Nu*strainf_exact(connu_e);
     
       % Energy norm solid displacement
       eENus_num = eENus_num + (esh-ese)^2 * Quad.w(ip,1) * Jdet;
       eENus_den = eENus_den + (ese)^2 * Quad.w(ip,1) * Jdet;

       % Energy norm pressure
       eENp_num = eENp_num + (gph-gpe)^2 * Quad.w(ip,1) * Jdet;
       eENp_den = eENp_den + (gpe)^2 * Quad.w(ip,1) * Jdet;

       % Energy norm fluid displacement
       eENuf_num = eENuf_num + (efh-efe)^2 * Quad.w(ip,1) * Jdet;
       eENuf_den = eENuf_den + (efe)^2 * Quad.w(ip,1) * Jdet;

       % H1 norm solid displacement
       eH1us_num = eH1us_num + ((esh-ese)^2 + (ush-use)^2) * Quad.w(ip,1) * Jdet;
       eH1us_den = eH1us_den + ((ese)^2 + (use)^2) * Quad.w(ip,1) * Jdet;

       % H1 norm pressure
       eH1p_num = eH1p_num + ((gph-gpe)^2 + (ph-pe)^2) * Quad.w(ip,1) * Jdet;
       eH1p_den = eH1p_den + ((gpe)^2 + (pe)^2) * Quad.w(ip,1) * Jdet;

       % H1 norm fluid displacement
       eH1uf_num = eH1uf_num + ((efh-efe)^2 + (ufh-ufe)^2) * Quad.w(ip,1) * Jdet;
       eH1uf_den = eH1uf_den + ((efe)^2 + (ufe)^2) * Quad.w(ip,1) * Jdet;
    end
end

% square root
ErrorComp.eL2u = sqrt(eL2us_num/eL2us_den);
ErrorComp.eL2p = sqrt(eL2p_num/eL2p_den);
ErrorComp.eL2uf = sqrt(eL2uf_num/eL2uf_den);

ErrorComp.eENu = sqrt(eENus_num/eENus_den);
ErrorComp.eENp = sqrt(eENp_num/eENp_den);    
ErrorComp.eENuf = sqrt(eENuf_num/eENuf_den);

ErrorComp.eH1u = sqrt(eH1us_num/eH1us_den);
ErrorComp.eH1p = sqrt(eH1p_num/eH1p_den); 
ErrorComp.eH1uf = sqrt(eH1uf_num/eH1uf_den);

end