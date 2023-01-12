function [m_L2, m_e] = ManufacturedSolution1D_check(d1, d2, d3, s1, s2, s3, e1, e2, e3, Mesh1, Mesh2, Mesh3, Material, Control, Quad)
%MANUFACTUREDSOLUTION_CHECK Calculates the convergence rates
%   [m_L2, m_e] = ManufacturedSolution_check(d1, d2, d3, s1, s2, s3, Mesh1,
%   Mesh2, Mesh3) calculates the rates of convergence of the L2 error norm
%   and the energy norm for using a manufactured solution
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d1, d2, d3:             Displacement vectors from 3 runs with increasing mesh
%                           refinement
%   s1, s2, s3:             Nodal stress data from 3 runs with increasing mesh
%                           refinement
%   e1, e2, e3:             Nodal strain data from 3 runs with increasing mesh
%                           refinement
%   Mesh1, Mesh2, Mesh3:    Mesh data structure from 3 runs with increasing mesh
%                           refinement
%   1 - coarsest mesh, 2 - medium mesh, 3 - finest mesh

% Manufactured solution
%   ux = x^5 -x^4
%   exx = 5x^4 -4x^3
%   sxx = E * exx
% ------------------------------------------------------------------------
% Adapted from https://github.com/GCMLab (Acknowledgements: Bruce Gee)
% ------------------------------------------------------------------------

plot_on = 1; % turn plots on/off - debugging tool

% material parameters
E = Material.E;
nu = Material.nu;
% quadrature data
nq = Control.nqU^Mesh1.nsd; % total number of integration points

%% Step 1 - Loop through each element and calculate the L2 and e-norms
eL2 = zeros(3,1);
eEN = zeros(3,1);
h = zeros(3,1);

% loop over meshes
for sim = 1:3
    switch sim
        case 1
            d = d1;
            s = s1;
            e = e1;
            Mesh = Mesh1;
        case 2
            d = d2;
            s = s2;
            e = e2;
            Mesh = Mesh2;
        case 3
            d = d3;
            s = s3;
            e = e3;
            Mesh = Mesh3;
    end
    
    % Mesh size
    h(sim) = max(Mesh.coords)/Mesh.ne;
    
    % Calculate exact solutions
    x = Mesh.coords;
    
    d_exact = @(x) x.^5 - x.^4;
    e_exact = 5*x.^4 - 4*x.^3;
    s_exact = E*e_exact;

    % Calculate error norms
    eL2_num = 0;
    eL2_den = 0;
    eEN_num = 0;
    eEN_den = 0;
    
    % loop over elements
    for i = 1:Mesh.ne
        conn_e = Mesh.conn(i,:);
        gcoords = Mesh.coords(conn_e,:);
        
        % loop over integration points
        for ip = 1:nq
            % N matrices
            N = getN(Mesh, Quad, ip);
            % N derivatives
            dN = getdN(Mesh, Quad, ip);
            % Jacobian matrix
            J = dN*gcoords;
            % Jacobian determinant
            Jdet = det(J);
            
            % approximated displacement at quadrature point
            uxh_p = N*d(Mesh.xdofs(conn_e));
            % exact displacement at quadrature point
%             uxe_p = N*d_exact(Mesh.xdofs(conn_e));
            uxe_p = d_exact(N*Mesh.coords(conn_e));
            % L2 norm
            eL2_num = eL2_num + (uxh_p - uxe_p) * (uxh_p - uxe_p) * Quad.w(ip,1) * Jdet;
            eL2_den = eL2_den + (uxe_p * uxe_p) * Quad.w(ip,1) * Jdet;
            
            % approximated stress and strain at quadrature point
            ehp = (e(conn_e,:).') * N';
            shp = (s(conn_e,:).') * N';
            % exact stress and strain
            eep = (e_exact(conn_e,:).') * N';
            sep = (s_exact(conn_e,:).') * N';
            % e norm
            eEN_num = eEN_num + (ehp-eep)'*(shp - sep) * Quad.w(ip,1) * Jdet;
            eEN_den = eEN_den + eep'*sep * Quad.w(ip,1) * Jdet;
        end
    end
    
    eL2(sim) = sqrt(eL2_num/eL2_den);
    eEN(sim) = sqrt(eEN_num/eEN_den);
    
end

% Determine slope of L2 norm
pL2 = polyfit(log(h), log(eL2),1);
m_L2 = pL2(1);

pEN = polyfit(log(h), log(eEN),1);
m_e = pEN(1);

%% Step 2 - Calculate the slope of each curve
if plot_on
    figure;
    loglog(h,eL2,'o', 'LineWidth', 1.5);
    hold on;
    loglog(h,exp(pL2(2))*h.^pL2(1), 'k', 'LineWidth', 1.5);
    hold off;
    xlabel('Mesh size');
    ylabel('L2-norm');
    title(sprintf('1D Convergence order of L2-norm: %.2f', m_L2));
    
    figure;
    loglog(h,eEN,'o', 'LineWidth', 1.5);
    hold on;
    loglog(h,exp(pEN(2))*h.^pEN(1), 'k', 'LineWidth', 1.5);
    hold off;
    xlabel('Mesh size');
    ylabel('e-norm');
    title(sprintf('1D Convergence order of e-norm: %.2f', m_e));
end

end