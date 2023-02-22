function [ErrorU] = ComputeMeshSizeErrorU(MeshU, Solution, Plot, Control)
% ------------------------------------------------------------------------
% Compute error norms
% L2-norm: uses displacement
% Energy norm: uses strain and stress
% H1-norm: uses displacement and strain
% ------------------------------------------------------------------------

%% Initialize variables
% number of elements
ne = MeshU.ne;

% mesh element size
ErrorU.h = max(MeshU.coords)/MeshU.ne;

% approximate solutions
du = Solution.u;
e = Solution.e;
s = Solution.s;

% exact solutions evaluated at nodes
du_exact = Plot.uan_space;
e_exact = Solution.e_an;
s_exact = Solution.s_an;

% exact solution with symbolic function
d_exact = Control.uan_symb;

% higher quadrature to evaluate error
Control.nqU = 16;
Quad = GlobalQuad(MeshU, Control);

% numerator and denominator terms
eL2u_num = 0;
eL2u_den = 0;

eENu_num = 0;
eENu_den = 0;

eH1u_num = 0;
eH1u_den = 0;

%% Compute error norm
% loop over elements
for i = 1:ne
    % element connectivity
    connu_e = MeshU.conn(i,:);
    % global coordinates
    gcoords_u = MeshU.coords(connu_e,:);
    
    % loop over integration points
    for ip = 1:Quad.nq
       % N matrices
       Nu = getN(MeshU, Quad, ip);
       % N derivatives
       dNu = getdN(MeshU, Quad, ip);
       % Jacobian matrix
       J = dNu*gcoords_u;
       % Jacobian determinant
       Jdet = det(J);
       
       % values at the quadrature point
       % approximate (FEM)
       uxh = Nu*du(MeshU.DOF(connu_e)');       
       % exact
       uxe = Nu*du_exact(MeshU.DOF(connu_e)');
       
       % with function handle
       uxe = eval(subs(d_exact, Nu*MeshU.coords(connu_e)));
       
       % L2 norm displacement
       eL2u_num = eL2u_num + (uxh - uxe)^2 * Quad.w(ip,1) * Jdet;
       eL2u_den = eL2u_den + (uxe)^2 * Quad.w(ip,1) * Jdet;
              
       % approximate (FEM)
       eh = Nu*e(connu_e);
       sh = Nu*s(connu_e);
       
       % exact
       ee = Nu*e_exact(connu_e);
       se = Nu*s_exact(connu_e);
       
%        ee = eval(subs(e_exact, Nu*MeshU.coords(connu_e)));
%        se = eval(subs(s_exact, Nu*MeshU.coords(connu_e)));
     
       % energy norm
       eENu_num = eENu_num + (eh-ee) * (sh-se) * Quad.w(ip,1) * Jdet;
       eENu_den = eENu_den + ee * se * Quad.w(ip,1) * Jdet;

       % H1 norm
       eH1u_num = eH1u_num + ((eh-ee)^2 + (uxh-uxe)^2) * Quad.w(ip,1) * Jdet;
       eH1u_den = eH1u_den + ((ee)^2 + (uxe)^2) * Quad.w(ip,1) * Jdet;
    end
end

% square root
ErrorU.eL2u = sqrt(eL2u_num/eL2u_den);
ErrorU.eENu = sqrt(eENu_num/eENu_den);
ErrorU.eH1u = sqrt(eH1u_num/eH1u_den);

end