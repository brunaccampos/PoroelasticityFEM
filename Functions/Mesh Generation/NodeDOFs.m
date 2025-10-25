% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function Mesh = NodeDOFs(Mesh)
% Updates the structure array containing mesh information with an array
% of degree of freedom indices
% Acknowledgemnents: Chris Ladubec, Matin Parchei-Esfahani

switch Mesh.field
    case 'u'
        Mesh.nDOFe = Mesh.nne*Mesh.nsd;         % number of DOF per element
        Mesh.nDOF = Mesh.nn*Mesh.nsd;           % total number of DOF
        Mesh.DOF = zeros(Mesh.nn, Mesh.nsd);

        for sd = 1:Mesh.nsd
            Mesh.DOF(:,sd) = (sd : Mesh.nsd : (Mesh.nDOF-(Mesh.nsd-sd)))';
        end

    case {'p', 'n'}
        Mesh.nDOFe = Mesh.nne;         % number of DOF per element
        Mesh.nDOF = Mesh.nn;           % total number of DOF
        Mesh.DOF = (1:Mesh.nDOF).';
end

end