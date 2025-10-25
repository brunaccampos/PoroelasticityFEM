% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function [L] = getGatherMatrix (rows,columns,e)
% 'Gather' matrix

L = zeros(rows, columns);
for i = 1:rows
    for j = 1:columns
        if j == (rows-1)*(e-1)+i
            L(i,j) = 1;
        end
    end
end

end