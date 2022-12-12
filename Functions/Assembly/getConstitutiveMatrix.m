function [C] = getConstitutiveMatrix(Material, Mesh)
% Compute constitutive matrix
nsd = Mesh.nsd;

switch nsd
    case 1
        C = Material.E;
    case 2
        C = zeros(3,3);
        aux = Material.E/(1-Material.nu^2);
        C(1,1) = 1;
        C(2,2) = 1;
        C(3,3) = (1-Material.nu)/2;
        C(1,2) = Material.nu;
        C(2,1) = Material.nu;
        C = C*aux;
end