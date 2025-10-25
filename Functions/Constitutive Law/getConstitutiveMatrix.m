% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [C] = getConstitutiveMatrix(nMat, Material, Mesh)
% Return linear elastic constitutive matrix for solid constituent

% Material properties
E = Material.M(nMat).E;
nu = Material.M(nMat).nu;

% Compute constitutive matrix
nsd = Mesh.nsd;

switch nsd
    case 1
    switch Material.constLaw
        case 'PlaneStress'
            C = E;
        case 'PlaneStrain'
            C = E*(1-nu)/((1+nu)*(1-2*nu));
    end
    case 2
        switch Material.constLaw
            case 'PlaneStress'
                C = zeros(3,3);
                aux = E/(1-nu^2);
                C(1,1) = 1;
                C(2,2) = 1;
                C(3,3) = (1-nu)/2;
                C(1,2) = nu;
                C(2,1) = nu;
                C = C*aux;
            case 'PlaneStrain'
                C = zeros(3,3);
                aux = E/((1+nu) * (1-2*nu));
                C(1,1) = 1 - nu;
                C(2,2) = 1 - nu;
                C(3,3) = (1-2*nu)/2;
                C(1,2) = nu;
                C(2,1) = nu;
                C = C*aux;
        end
end