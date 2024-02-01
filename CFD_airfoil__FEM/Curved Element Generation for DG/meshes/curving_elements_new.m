function [node_loc, node_add, elem_rem] = curving_elements_new(in_gri, Q)
%curving_elements takes a gri file and curves the elements to match the
% geometry of the airfoil
%INPUT: in_gri is the input gri file to be edited
%       Q is the value to which the curving should occur
%OUTPUT: node_loc array of the new node locations for the curved elements
%        node_add array of the new nodes for each element along the
%        boundary


    [B2E, E2N, V] = Task2_Task3(in_gri);
    
    elem_rem = zeros(length(B2E), 1);
    start_idx = length(V);
    
    node_loc = zeros(10*length(B2E), 2);  % overallocate
    ff = length(find(B2E(:,3) == -1));
    
    if Q==2
        node_add = zeros(length(B2E)-ff, 6);
    end
    if Q==3
        node_add = zeros(length(B2E)-ff, 10);
    end

    j = 1;      % counting index for adding the new node locations
    
    for i=1:length(E2N)
         
        non = find(B2E(:,1) == i);
        
        if non ~= 0
            continue
            
        else
            nodes = E2N(i, :);
            str_nds_1 = straight_nodes(V(nodes(1),:), V(nodes(2),:), Q);
            str_nds_2 = straight_nodes(V(nodes(2),:), V(nodes(3),:), Q);
            str_nds_3 = straight_nodes(V(nodes(3),:), V(nodes(1),:), Q);
            
            if Q == 2
                node_loc(j, :) = str_nds_1;
                node_loc(j+1, :) = str_nds_3;
                node_loc(j+2, :) = str_nds_2;
                
                node_add(i,:) = [nodes(1), start_idx+1, nodes(2), start_idx+2, start_idx+3, nodes(3)];
            end
            
            if Q == 3
                center_node_x = mean([str_nds_1(:,1), str_nds_2(:,1), str_nds_3(:,1)], 'all');
                center_node_y = mean([str_nds_1(:,2), str_nds_2(:,2), str_nds_3(:,2)], 'all');

                node_loc(j,:) = str_nds_1(1,:);
                node_loc(j+1,:) = str_nds_1(2,:);
                node_loc(j+2,:) = str_nds_3(2,:);
                node_loc(j+3,1) = center_node_x;
                node_loc(j+3,2) = center_node_y;
                node_loc(j+4,:) = str_nds_2(1,:);
                node_loc(j+5,:) = str_nds_3(1,:);
                node_loc(j+6,:) = str_nds_2(2,:);

                node_add(i,:) = [nodes(1), start_idx+1, start_idx+2, nodes(2), start_idx+3, start_idx+4, start_idx+5, start_idx+6, start_idx+7, nodes(3)];
            end    
            
            if Q == 2
                start_idx = start_idx+3;
                j = j+3;
                
            end
            if Q==3
                start_idx = start_idx+7;
                j = j+7;
                
            end
                      
        end
        
    end
        
    
    for i=1:length(B2E)
    
        if B2E(i, 3) ~= -1
            
            if B2E(i,3) == -2
                flag = 3;
            elseif B2E(i,3) == -3
                flag = 2;
            else
                flag = 1;
            end
            
            elem = B2E(i, 1);
            nodes = E2N(elem, :);
            
            elem_rem(i) = elem;
            
            curve_num = B2E(i,2);
            if curve_num == 1
                curve = curve_nodes(V(nodes(2),:), V(nodes(3),:), Q, flag);
                straight_nodes_1 = straight_nodes(V(nodes(3),:), V(nodes(1),:), Q);
                straight_nodes_2 = straight_nodes(V(nodes(1),:), V(nodes(2),:), Q);
                
                if Q == 2
                    node_loc(j+1,:) = straight_nodes_2;
                    node_loc(j+2,:) = straight_nodes_1;
                    node_loc(j,:) = curve;

                    node_add(elem, :) = [nodes(2), start_idx+1, nodes(3), start_idx+2, start_idx+3, nodes(1)];
                end
                
                if Q == 3
                    center_node_x = mean([curve(:,1), straight_nodes_1(:,1), straight_nodes_2(:,1)], 'all');
                    center_node_y = mean([curve(:,2), straight_nodes_1(:,2), straight_nodes_2(:,2)], 'all');
                
                    node_loc(j,:) = curve(1,:);
                    node_loc(j+1,:) = curve(2,:);
                    node_loc(j+2,:) = straight_nodes_2(2,:);
                    node_loc(j+3,1) = center_node_x;
                    node_loc(j+3,2) = center_node_y;
                    node_loc(j+4,:) = straight_nodes_2(1,:);
                    node_loc(j+5,:) = straight_nodes_1(1,:);
                    node_loc(j+6,:) = straight_nodes_1(2,:);
                    
                    node_add(elem,:) = [nodes(2), start_idx+1, start_idx+2, nodes(3), start_idx+3, start_idx+4, start_idx+5, start_idx+6, start_idx+7, nodes(1)];
                end
                
            elseif curve_num == 2
                curve = curve_nodes(V(nodes(3),:), V(nodes(1),:), Q, flag);
                straight_nodes_1 = straight_nodes(V(nodes(1),:), V(nodes(2),:), Q);
                straight_nodes_2 = straight_nodes(V(nodes(2),:), V(nodes(3),:), Q);
                
                if Q == 2
                    node_loc(j+1,:) = straight_nodes_2;
                    node_loc(j+2,:) = straight_nodes_1;
                    node_loc(j,:) = curve;

                    node_add(elem, :) = [nodes(3), start_idx+1, nodes(1), start_idx+2, start_idx+3, nodes(2)];
                end
                
                if Q == 3
                    center_node_x = mean([curve(:,1), straight_nodes_1(:,1), straight_nodes_2(:,1)], 'all');
                    center_node_y = mean([curve(:,2), straight_nodes_1(:,2), straight_nodes_2(:,2)], 'all');
                
                    node_loc(j,:) = curve(1,:);
                    node_loc(j+1,:) = curve(2,:);
                    node_loc(j+2,:) = straight_nodes_2(2,:);
                    node_loc(j+3,1) = center_node_x;
                    node_loc(j+3,2) = center_node_y;
                    node_loc(j+4,:) = straight_nodes_2(1,:);
                    node_loc(j+5,:) = straight_nodes_1(1,:);
                    node_loc(j+6,:) = straight_nodes_1(2,:);
                    
                    node_add(elem,:) = [nodes(3), start_idx+1, start_idx+2, nodes(1), start_idx+3, start_idx+4, start_idx+5, start_idx+6, start_idx+7, nodes(2)];
                end
                
            else
                curve = curve_nodes(V(nodes(1),:), V(nodes(2),:), Q, flag);
                straight_nodes_1 = straight_nodes(V(nodes(2),:), V(nodes(3),:), Q);
                straight_nodes_2 = straight_nodes(V(nodes(3),:), V(nodes(1),:), Q);
                
                if Q == 2
                    node_loc(j+1,:) = straight_nodes_2;
                    node_loc(j+2,:) = straight_nodes_1;
                    node_loc(j,:) = curve;

                    node_add(elem, :) = [nodes(1), start_idx+1, nodes(2), start_idx+2, start_idx+3, nodes(3)];
                end
                
                if Q == 3
                    center_node_x = mean([curve(:,1), straight_nodes_1(:,1), straight_nodes_2(:,1)], 'all');
                    center_node_y = mean([curve(:,2), straight_nodes_1(:,2), straight_nodes_2(:,2)], 'all');
                
                    node_loc(j,:) = curve(1,:);
                    node_loc(j+1,:) = curve(2,:);
                    node_loc(j+2,:) = straight_nodes_2(2,:);
                    node_loc(j+3,1) = center_node_x;
                    node_loc(j+3,2) = center_node_y;
                    node_loc(j+4,:) = straight_nodes_2(1,:);
                    node_loc(j+5,:) = straight_nodes_1(1,:);
                    node_loc(j+6,:) = straight_nodes_1(2,:);
                    
                    node_add(elem,:) = [nodes(1), start_idx+1, start_idx+2, nodes(2), start_idx+3, start_idx+4, start_idx+5, start_idx+6, start_idx+7, nodes(3)];
                end
                
            end
            
            if Q == 2
                start_idx = start_idx+3;
                j = j+3;
                
            end
            if Q==3
                start_idx = start_idx+7;
                j = j+7;
                
            end
            
        
        else
            elem = B2E(i, 1);
            nodes = E2N(elem, :);
            
            str_nds_1 = straight_nodes(V(nodes(1),:), V(nodes(2),:), Q);
            str_nds_2 = straight_nodes(V(nodes(2),:), V(nodes(3),:), Q);
            str_nds_3 = straight_nodes(V(nodes(3),:), V(nodes(1),:), Q);
            
            if Q == 2
                node_loc(j, :) = str_nds_1;
                node_loc(j+1, :) = str_nds_3;
                node_loc(j+2, :) = str_nds_2;
                
                node_add(elem,:) = [nodes(1), start_idx+1, nodes(2), start_idx+2, start_idx+3, nodes(3)];
            end
            
            if Q == 3
                center_node_x = mean([str_nds_1(:,1), str_nds_2(:,1), str_nds_3(:,1)], 'all');
                center_node_y = mean([str_nds_1(:,2), str_nds_2(:,2), str_nds_3(:,2)], 'all');

                node_loc(j,:) = str_nds_1(1,:);
                node_loc(j+1,:) = str_nds_1(2,:);
                node_loc(j+2,:) = str_nds_3(2,:);
                node_loc(j+3,1) = center_node_x;
                node_loc(j+3,2) = center_node_y;
                node_loc(j+4,:) = str_nds_2(1,:);
                node_loc(j+5,:) = str_nds_3(1,:);
                node_loc(j+6,:) = str_nds_2(2,:);

                node_add(elem,:) = [nodes(1), start_idx+1, start_idx+2, nodes(2), start_idx+3, start_idx+4, start_idx+5, start_idx+6, start_idx+7, nodes(3)];
            end    
            
            if Q == 2
                start_idx = start_idx+3;
                j = j+3;
                
            end
            if Q==3
                start_idx = start_idx+7;
                j = j+7;
                
            end
            
        end
    end
    
    node_loc = node_loc(1:j-1,:);
    elem_rem = nonzeros(elem_rem);

end

