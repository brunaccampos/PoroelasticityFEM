function [ErrorComp] = ComputeMeshSizeError(MeshU, MeshP, Solution, Plot, Control)
% ------------------------------------------------------------------------
% Compute displacement and pressure errors related to the mesh size
% ------------------------------------------------------------------------

%% Initialize variables
% number of elements
ne = MeshU.ne;

% mesh element size
switch MeshP.nsd
    case 1
        ErrorComp.h = max(MeshP.coords)/MeshP.ne;
        
        MeshU.xdofs = MeshU.DOF;
        MeshU.ydofs = zeros(length(MeshU.DOF),1);
        MeshP.xdofs = MeshP.DOF;
        MeshP.ydofs = zeros(length(MeshP.DOF),1);
    case 2
        gcoordsp = MeshP.coords(MeshP.conn(1,:),:);
        ErrorComp.h = sqrt(polyarea(gcoordsp(:,1),gcoordsp(:,2)));
end

% approximate solutions
du = Solution.u;
dp = Solution.p;
e = Solution.e;
s = Solution.s;
q = Solution.q;

% exact solutions evaluated at nodes
du_exact = Plot.uan_space;
dp_exact = Plot.pan_space;
e_exact = Solution.e_an;
s_exact = Solution.s_an;
q_exact = Solution.q_an;

% d_exact = @(x) x.^5 - x.^4;
% e_exact = @(x) 5*x.^4 - 4*x.^3;
% s_exact = e_exact;

% p_exact = @(x) x.^5 - x.^4;
% qex_exact = @(x)-( 5*x.^4 - 4*x.^3);

% syms x
% aux=0;
% % x = MeshP.coords;
% N=1000;
% % for k=1:N
% %     aux = aux + (1/k)*exp(-3*k^2*pi()^2*10)*sin(k*pi()*x);
% % end
% % 
% % p_exact = 100 - 50*x - (100/pi()) * aux;
% p_exact = 100 - 50*x;
% 
        
        
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
       Np = getN(MeshP, Quad, ip);
       Nu = getN(MeshU, Quad, ip);
       % N derivatives
       dNp = getdN(MeshP, Quad, ip);
       % Jacobian matrix
       J = dNp*gcoords_p;
       % Jacobian determinant
       Jdet = det(J);
       
       % values at the quadrature point
       switch MeshP.nsd
           case 1 % 1D
               % approximate (FEM)
               uxh = Nu*du(MeshU.xdofs(connu_e)');
               uyh = zeros(length(uxh));
               
               pxh = Np*dp(MeshP.xdofs(connp_e)');
               pyh = zeros(length(pxh));
               
               % exact
               uxe = Nu*du_exact(MeshU.xdofs(connu_e)');
               uye = zeros(length(uxe));
               
               pxe = Np*dp_exact(MeshP.xdofs(connp_e)');
               pye = zeros(length(pxe));
               
               % with function handle
%                uxe = eval(subs(d_exact, Nu*MeshU.coords(connu_e)));
%                pxe = eval(subs(p_exact, Np*MeshP.coords(connp_e)));
               
           case 2 % 2D
               % approximate (FEM)
               uxh = Nu*du(MeshU.xdofs(connu_e)');
               uyh = Nu*du(MeshU.ydofs(connu_e)');
               
               pxh = Np*dp(MeshP.xdofs(connp_e)');
               pyh = Np*dp(MeshP.ydofs(connp_e)');
               
               % exact
               uxe = Nu*du_exact(MeshU.xdofs(connu_e)');
               uye = Nu*du_exact(MeshU.ydofs(connu_e)');
               
               pxe = Np*dp_exact(MeshP.xdofs(connp_e)');
               pye = Np*dp_exact(MeshP.ydofs(connp_e)'); 
       end
       
       % L2 norm displacement
       eL2u_num = eL2u_num + [uxh - uxe, uyh - uye] * [uxh - uxe; uyh - uye] * Quad.w(ip,1) * Jdet;
       eL2u_den = eL2u_den + [uxe, uye] * [uxe; uye] * Quad.w(ip,1) * Jdet;
       
       % L2 norm pressure
       eL2p_num = eL2p_num + [pxh - pxe, pyh - pye] * [pxh - pxe; pyh - pye] * Quad.w(ip,1) * Jdet;
       eL2p_den = eL2p_den + [pxe, pye] * [pxe; pye] * Quad.w(ip,1) * Jdet;
       
       % approximate (FEM)
       eh = Nu*e(connu_e);
       sh = Nu*s(connu_e);
       qh = Np*q(connp_e);
       
       % exact
       ee = Nu*e_exact(connu_e);
       se = Nu*s_exact(connu_e);
       qe = Np*q_exact(connp_e);
       
%        ee = eval(subs(e_exact, Nu*MeshU.coords(connu_e)));
%        se = eval(subs(s_exact, Nu*MeshU.coords(connu_e)));

%        qe = eval(subs(qex_exact, Np*MeshP.coords(connp_e)));
     
       % energy norm
       eENu_num = eENu_num + (eh-ee) * (sh-se) * Quad.w(ip,1) * Jdet;
       eENu_den = eENu_den + ee * se * Quad.w(ip,1) * Jdet;
       eENp_num = eENp_num + (qh-qe) * (qh-qe) * Quad.w(ip,1) * Jdet;
       eENp_den = eENp_den + qe * qe * Quad.w(ip,1) * Jdet;


       % H1 norm
       eH1u_num = eH1u_num + ((eh-ee)^2 + (uxh - uxe)^2) * Quad.w(ip,1) * Jdet;
       eH1u_den = eH1u_den + ((ee)^2 + (uxe)^2) * Quad.w(ip,1) * Jdet;
       eH1p_num = eH1p_num + ((qh-qe)^2 + (pxh - pxe)^2) * Quad.w(ip,1) * Jdet;
       eH1p_den = eH1p_den + ((qe)^2 + (pxe)^2) * Quad.w(ip,1) * Jdet;
    end
end

% square root
ErrorComp.eL2u = sqrt(eL2u_num/eL2u_den);
ErrorComp.eL2p = sqrt(eL2p_num/eL2p_den);

ErrorComp.eENp = sqrt(eENp_num/eENp_den);    
ErrorComp.eENu = sqrt(eENu_num/eENu_den);

ErrorComp.eH1p = sqrt(eH1p_num/eH1p_den);    
ErrorComp.eH1u = sqrt(eH1u_num/eH1u_den);

end