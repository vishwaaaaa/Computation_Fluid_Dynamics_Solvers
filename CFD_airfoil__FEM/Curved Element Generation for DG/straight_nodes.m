function [new_locs] = straight_nodes(node1_loc, node2_loc, Q)
%straight_nodes determines the new locations of the additional nodes in
%each element, using counter clockwise ordering
% INPUT : node1_loc = [1, 2] array of the location of the first node
%         node2_loc = [1, 2] array of the location of the second node
%         Q the amount of curving, only supports Q=2 or 3
% OUTPUT: new_locs = [Q-1, 2] array of the additional node locations

new_locs = zeros(Q-1,2);

if Q == 2,
   
    new_locs(1, 1) = (node1_loc(1) + node2_loc(1))/2;
    new_locs(1, 2) = (node1_loc(2) + node2_loc(2))/2;
    
end

if Q == 3,
    
    new_locs(1, 1) = (node1_loc(1) + node2_loc(1))/3;
    new_locs(1, 2) = (node1_loc(2) + node2_loc(2))/3;
    
    new_locs(2, 1) = 2*(node1_loc(1) + node2_loc(1))/3;
    new_locs(2, 2) = 2*(node1_loc(2) + node2_loc(2))/3;
    
end

end

