% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Bvoigt] = getBVoigt(Mesh, B)
% return the B matrix in Voigt form
% Acknowledgements: Matin Parchei Esfahani

nsd = Mesh.nsd;

switch nsd
    case 1
        Bvoigt = B;
    case 2
        if Mesh.field == 'p' || Mesh.field == 'n'
            Bvoigt = B;
        else
            Bvoigt = zeros(3, Mesh.nsd*Mesh.nne);

            Bvoigt(1,1:2:end) = B(1,:);
            Bvoigt(2,2:2:end) = B(2,:);
            Bvoigt(3,1:2:end) = B(2,:);
            Bvoigt(3,2:2:end) = B(1,:);
        end
end

end