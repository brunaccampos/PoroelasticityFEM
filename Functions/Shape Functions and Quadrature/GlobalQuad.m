function [Quad] = GlobalQuad(Mesh, Control)
% Select points of global quadrature

% number of spatial dimensions
nsd = Mesh.nsd;
% quadrature order
nq = Control.nq;

% quadrature coordinates
Quad.p = zeros(nq^nsd,nsd);
% quadratrure weights
Quad.w = zeros(nq^nsd,1);

pt = zeros(nq,1);
wt = zeros(nq,1);

switch nq
    case 2
        % coords
        pt(1,1) = -0.5773502692;
        pt(2,1) = 0.5773502692;
        % weights
        wt(:,1) = 1;
end

n = 1; % counter

% store points according to the spatial dimension
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