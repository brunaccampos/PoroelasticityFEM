function [Quad] = GlobalQuad(Control)

% quadrature coordinates
Quad.csi = zeros(Control.nq,1);

% quadratrure weights
Quad.w_csi = zeros(Control.nq,1);

switch Control.nq
    case 2
        Quad.csi(1,1) = -0.5773502692;
        Quad.csi(2,1) = 0.5773502692;

        Quad.w_csi(:,1) = 1;
end 

end