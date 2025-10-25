% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [d] = MatrixInvert(K,f,parallel_number)
% Solves the linear system d = inv(K)*f

if parallel_number == 1
    d = K\f;
else
    % Check if < 200k dofs
        if length(f) < 2e5
            fprintf(' Careful! Matrix is being inverted in parallel at < 2e5 degrees of freedom. \n At < 2e5 degrees of freedom, single core is generally more efficient.\n')
        end
    % Create parallel processing pool
        if isempty(gcp('nocreate'))
            parpool(parallel_number);
        end
    % distribute to pool
        K = distributed(K);
        f = distributed(f);
    % Invert matrix
        d = K\f;
    % Return distributed to double
        d = gather(d);
end

