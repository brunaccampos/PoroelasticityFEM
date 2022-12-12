function [Quad] = GlobalQuad(Mesh, Control)
% Select points of global quadrature
%   --------------------------------------------------------------------
%   References
%   --------------------------------------------------------------------
%   [1] Fish, J., & Belytschko, T. (2007). A first course in finite
%       elements. Chichester; Hoboken, NJ: John Wiley & Sons.

%   Adapted from: https://github.com/GCMLab (Acknowledgements: Jack Chessa)
%   --------------------------------------------------------------------

% number of spatial dimensions
nsd = Mesh.nsd;

% quadrature order
switch Mesh.field
    case 'u'
        nq = Control.nqU;
    otherwise
        nq = Control.nqP;
end

switch Mesh.type
    case {'T3', 'T6'}
        %% 2D triangular elements
        switch nq
            case 1
                % coords
                Quad.p = [1/3, 1/3];
                % weight
                Quad.w = 1/2;
            case 2
                Quad.p = zeros(3,2);
                Quad.w = zeros(3,1);
                % coords
                Quad.p(1,:) = [1/6, 1/6];
                Quad.p(2,:) = [4/6, 1/6];
                Quad.p(3,:) = [1/6, 4/6];
                % weights
                Quad.w(:) = 1/6;
        end
    otherwise
        %% 1D and 2D quadrilateral elements
        % quadrature coordinates
        Quad.p = zeros(nq^nsd,nsd);
        % quadratrure weights
        Quad.w = zeros(nq^nsd,1);
        
        pt = zeros(nq,1);
        wt = zeros(nq,1);
        switch nq
            case 1
                % coords
                pt(1) = 0;
                % weights
                wt(1) = 2;
            case 2
                % coords
                pt(1) = 0.577350269189626;
                pt(2) = -0.577350269189626;
                % weights
                wt(:) = 1;
            case 3
                % coords
                pt(1) =  0.774596669241483;
                pt(2) = -0.774596669241483;
                pt(3) =  0.000000000000000;
                % weights
                wt(1) = 0.555555555555556;
                wt(2) = 0.555555555555556;
                wt(3) = 0.888888888888889;
            case 4
                % coords
                pt(1) =  0.861134311594053;
                pt(2) = -0.861134311594053;
                pt(3) =  0.339981043584856;
                pt(4) = -0.339981043584856;
                % weights
                wt(1) = 0.347854845137454;
                wt(2) = 0.347854845137454;
                wt(3) = 0.652145154862546;
                wt(4) = 0.652145154862546;
            case 5
                % coords
                pt(1) =  0.906179845938664;
                pt(2) = -0.906179845938664;
                pt(3) =  0.538469310105683;
                pt(4) = -0.538469310105683;
                pt(5) =  0.000000000000000;
                % weights
                wt(1) = 0.236926885056189;
                wt(2) = 0.236926885056189;
                wt(3) = 0.478628670499366;
                wt(4) = 0.478628670499366;
                wt(5) = 0.568888888888889;
            case 6
                % coords
                pt(1) =  0.932469514203152;
                pt(2) = -0.932469514203152;
                pt(3) =  0.661209386466265;
                pt(4) = -0.661209386466265;
                pt(5) =  0.238619186003152;
                pt(6) = -0.238619186003152;
                % weights
                wt(1) = 0.171324492379170;
                wt(2) = 0.171324492379170;
                wt(3) = 0.360761573048139;
                wt(4) = 0.360761573048139;
                wt(5) = 0.467913934572691;
                wt(6) = 0.467913934572691;
            case 7
                % coords
                pt(1) =  0.949107912342759;
                pt(2) = -0.949107912342759;
                pt(3) =  0.741531185599394;
                pt(4) = -0.741531185599394;
                pt(5) =  0.405845151377397;
                pt(6) = -0.405845151377397;
                pt(7) =  0.000000000000000;
                % weights
                wt(1) = 0.129484966168870;
                wt(2) = 0.129484966168870;
                wt(3) = 0.279705391489277;
                wt(4) = 0.279705391489277;
                wt(5) = 0.381830050505119;
                wt(6) = 0.381830050505119;
                wt(7) = 0.417959183673469;
            case 8
                % coords
                pt(1) =  0.960289856497536;
                pt(2) = -0.960289856497536;
                pt(3) =  0.796666477413627;
                pt(4) = -0.796666477413627;
                pt(5) =  0.525532409916329;
                pt(6) = -0.525532409916329;
                pt(7) =  0.183434642495650;
                pt(8) = -0.183434642495650;
                % weights
                wt(1) = 0.101228536290376;
                wt(2) = 0.101228536290376;
                wt(3) = 0.222381034453374;
                wt(4) = 0.222381034453374;
                wt(5) = 0.313706645877887;
                wt(6) = 0.313706645877887;
                wt(7) = 0.362683783378362;
                wt(8) = 0.362683783378362;
        end
end
n = 1; % counter

% store points according to the spatial dimension
if strcmp(Mesh.type,'T3') == 0 && strcmp(Mesh.type,'T6') == 0
    if(nsd == 1)
        Quad.p = pt;
        Quad.w = wt;
    elseif(nsd == 2)
        for i = 1:nq
            for j = 1:nq
                Quad.p(n,:) = [pt(i,1), pt(j,1)];
                Quad.w(n) = wt(i,1)*wt(j,1);
                n = n+1;
            end
        end
    end
end

% number of quadrature points
Quad.nq = length(Quad.w);