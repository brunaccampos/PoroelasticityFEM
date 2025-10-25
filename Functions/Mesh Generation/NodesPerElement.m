% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function nne = NodesPerElement(etype)
% Number of nodes per element type

switch etype
    case 'L2'
        nne = 2;
    case 'L3'
        nne = 3;
    case 'L4'
        nne = 4;
    case 'T3'
        nne = 3;
    case 'T6'
        nne = 6;
    case 'Q4'
        nne = [2;2];
    case 'Q9'
        nne = [3;3];
end

end

