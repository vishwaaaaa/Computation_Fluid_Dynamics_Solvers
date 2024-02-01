function [new_locs] = curve_nodes(node1_loc, node2_loc, Q, flag)
%curve_nodes snaps the new nodes for the curved edges of the mesh to lie on
%the airfoil, using counter clockwise counting
%INPUT: node1_loc = [1,2] array of the first node location on curved edge
%       node2_loc = [1,2] array of the second node location on curved edge
%       Q is the curving order, only supports 2 and 3
%       flag is the element on which the curved edge is on
%OUTPUT: new_locs array of the new node locations

new_locs = zeros(Q-1,2);

if Q == 2
    
    straight_node_loc = (node1_loc + node2_loc)/2;
    V_mod = nspline(straight_node_loc, flag); 
    new_locs = V_mod;
    
end

if Q == 3
    
    str_loc_1 = [(node1_loc(1,1) + node2_loc(1,1))/3, (node1_loc(1,2) + node2_loc(1,2))/3];
    V_mod_1 = nspline(str_loc_1, flag);
    new_locs(1,:) = V_mod_1;
    
    str_loc_2 = [2*(node1_loc(1,1) + node2_loc(1,1))/3, 2*(node1_loc(1,2) + node2_loc(1,2))/3];
    V_mod_2 = nspline(str_loc_2, flag);
    new_locs(2,:) = V_mod_2;

end    

end

