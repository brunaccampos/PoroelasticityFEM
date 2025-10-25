% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [Nvoigt] = getNVoigt(Mesh, N)
% Return the shape functions matrix in Voigt form
% Acknowledgements: Matin Parchei Esfahani

% number of spatial dimensions
nsd = Mesh.nsd;

switch nsd
    case 1
        Nvoigt = N;
    case 2
        if Mesh.field == 'p' || Mesh.field == 'n'
            Nvoigt = N;
        else
            n = size(N,2);
            I = eye(Mesh.nsd);
            Nvoigt = [];
            for i = 1:n
                Nvoigt = [Nvoigt, I*N(i)];
            end
        end
end

end