function [ErrorComp] = ComputeMeshSizeError(MeshU, MeshP, Solution, Plot)
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

% exact solutions evaluated at nodes
du_exact = Plot.uan_space;
dp_exact = Plot.pan_space;

% higher quadrature to evaluate error
Control.nqU = 8;
Quad = GlobalQuad(MeshU, Control);

% numerator and denominator terms
eL2u_num = 0;
eL2u_den = 0;
eL2p_num = 0;
eL2p_den = 0;

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
               uxh_u = Nu*du(MeshU.xdofs(connu_e)');
               uyh_u = zeros(length(uxh_u));
               
               uxh_p = Np*dp(MeshP.xdofs(connp_e)');
               uyh_p = zeros(length(uxh_p));
               
               % exact
               uxe_u = Nu*du_exact(MeshU.xdofs(connu_e)');
               uye_u = zeros(length(uxe_u));
               
               uxe_p = Np*dp_exact(MeshP.xdofs(connp_e)');
               uye_p = zeros(length(uxe_p));
               
           case 2 % 2D
               % approximate (FEM)
               uxh_u = Nu*du(MeshU.xdofs(connu_e)');
               uyh_u = Nu*du(MeshU.ydofs(connu_e)');
               
               uxh_p = Np*dp(MeshP.xdofs(connp_e)');
               uyh_p = Np*dp(MeshP.ydofs(connp_e)');
               
               % exact
               uxe_u = Nu*du_exact(MeshU.xdofs(connu_e)');
               uye_u = Nu*du_exact(MeshU.ydofs(connu_e)');
               
               uxe_p = Np*dp_exact(MeshP.xdofs(connp_e)');
               uye_p = Np*dp_exact(MeshP.ydofs(connp_e)'); 
       end
       
       % L2 norm displacement
       eL2u_num = eL2u_num + [uxh_u - uxe_u, uyh_u - uye_u] * [uxh_u - uxe_u; uyh_u - uye_u] * Quad.w(ip,1) * Jdet;
       eL2u_den = eL2u_den + [uxe_u, uye_u] * [uxe_u; uye_u] * Quad.w(ip,1) * Jdet;
       
       % L2 norm pressure
       eL2p_num = eL2p_num + [uxh_p - uxe_p, uyh_p - uye_p] * [uxh_p - uxe_p; uyh_p - uye_p] * Quad.w(ip,1) * Jdet;
       eL2p_den = eL2p_den + [uxe_p, uye_p] * [uxe_p; uye_p] * Quad.w(ip,1) * Jdet;
    end
end

% square root
ErrorComp.eL2u = sqrt(eL2u_num/eL2u_den);
ErrorComp.eL2p = sqrt(eL2p_num/eL2p_den);

end