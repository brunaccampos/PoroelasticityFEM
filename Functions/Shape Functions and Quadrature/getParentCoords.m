% SPDX-FileCopyrightText: Copyright (c) 2022-2024 Bruna Campos
% SPDX-License-Identifier: GPL-3.0-or-later

function xi_node = getParentCoords(node_num, type)
% Parent coordinates of the nodal position specified

switch type
    case 'L2'
        if node_num == 1
            xi_node = -1;
        elseif node_num == 2
            xi_node = 1; 
        end
    case 'L3'
        if node_num == 1
            xi_node = -1;
        elseif node_num == 2
            xi_node = 0;
        elseif node_num == 3
            xi_node = 1;
        end
    case 'L4'
        if node_num == 1
            xi_node = -1;
        elseif node_num == 2
            xi_node = -1/3;
        elseif node_num == 3
            xi_node = 1/3;
        elseif node_num == 4
            xi_node = 1;
        end
    case 'T3'
        if node_num == 1
            xi_node = [0,0];
        elseif node_num == 2
            xi_node = [1,0];
        elseif node_num == 3
            xi_node = [0,1];
        end
    case 'T6'
        if node_num == 1
            xi_node = [0,0];
        elseif node_num == 2
            xi_node = [1,0];
        elseif node_num == 3
            xi_node = [0,1];
        elseif node_num == 4
            xi_node = [0.5,0];
        elseif node_num == 5
            xi_node = [0,0.5];
        elseif node_num == 6
            xi_node = [0.5,0.5];
        end
    case 'Q4'
        if node_num == 1
            xi_node = [-1 -1];
        elseif node_num == 2
            xi_node = [1 -1];
        elseif node_num == 3
            xi_node = [1 1];
        elseif node_num == 4
            xi_node = [-1 1];
        end
    case 'Q9'
        if node_num == 1
            xi_node = [-1 -1];
        elseif node_num == 2
            xi_node = [1 -1];
        elseif node_num == 3
            xi_node = [1 1];
        elseif node_num == 4
            xi_node = [-1 1];
        elseif node_num == 5
            xi_node = [0 -1];
        elseif node_num == 6
            xi_node = [1 0];
        elseif node_num == 7
            xi_node = [0 1];
        elseif node_num == 8
            xi_node = [-1 0];
        elseif node_num == 9
            xi_node = [0 0];
        end
    case 'B8'
        if node_num == 1
            xi_node = [-1 -1 -1];
        elseif node_num == 2
            xi_node = [1 -1 -1];
        elseif node_num == 3
            xi_node = [1 1 -1];
        elseif node_num == 4
            xi_node = [-1 1 -1];
        elseif node_num == 5
            xi_node = [-1 -1 1];
        elseif node_num == 6
            xi_node = [1 -1 1];
        elseif node_num == 7
            xi_node = [1 1 1];
        elseif node_num == 8
            xi_node = [-1 1 1];
        end
end

end

