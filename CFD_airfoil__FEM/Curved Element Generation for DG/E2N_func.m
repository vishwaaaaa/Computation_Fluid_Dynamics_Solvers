function [V, E2N, edge_info] = E2N_func(grifile)
% This function takes the standard .gri file for the mesh configuration and
% transfers it to the E2N format, saving the node locations and the
% boundary edge information
%
% INPUT : grifile str standard .gri file
% OUTPUT : E2N = [nelem x 3] array mapping triangles to nodes
%          V = [nnodes x 2] array of the location of each node
%          edge_info = [nB x 4] array with the information for each
%          boundary

fid = fopen(grifile, 'r');

% Read in nodes and node locations
A = fscanf(fid,'%d', 3);
nnode    = A(1);
nelemtot = A(2);
dim      = A(3);
V = zeros(nnode, dim);  % node locations
for inode = 1:nnode
  A = fscanf(fid, '%lf', 2);
  V(inode,:) = A(1:2)';
end

% Read through boundary info
A = fscanf(fid, '%d', 1);
nbfgrp = A(1);
for ibfgrp = 1:nbfgrp,
  fgets(fid);
  sline = fgets(fid);
  [nbface, nnode, title] = strread(sline, '%d %d %s');
  A = zeros(nbface);
  for ibface = 1:nbface,
    Ash = fscanf(fid, '%d', nnode);
    A(ibface) = Ash(1);
  end
  l{ibfgrp}=A;
end

% Doing the E2N calculation
curtot = 0;
E2N = zeros(nelemtot, 3);
while (curtot ~= nelemtot)
  fgets(fid);
  sline = fgets(fid);
  [nelem, p, sbasis] = strread(sline, '%d %d %s');
  switch sbasis{1}
    case 'TriLagrange'
      nnode = (p+1)*(p+2)/2;
      nedge = 3;
      fvec = zeros(3,p+1);
      fvec(1,:) = [1:(p+1)];
      v = p+1; d = p;
      for k=1:(p+1) 
          fvec(2,k) = v;
          v = v+d; 
          d = d-1; 
      end
      v = 1; 
      d = p+1;
      for k=1:(p+1) 
          fvec(3,k) = v; 
          v = v+d; 
          d = d-1; 
      end
    otherwise
      error('element type not understood');
  end
  
  E2N_int = zeros(nnode, nelem);
  for elem = 1:nelem
     E2N(elem, :) = fscanf(fid, '%d', nnode);  % HO nodes included
  end
  curtot = curtot + nelem;
end
fclose(fid);
nelem = nelemtot;

% Getting the information about the boundary
n_edge_inc = 3*nelem;
edge_info = zeros(n_edge_inc, 5);
for ielem=1:nelem
    edge_info((3*ielem)-2,1) = E2N(ielem,1);    %1st node of 1st edge - 1st node reported for the element
    edge_info((3*ielem)-2,2) = E2N(ielem,2);    %2nd node of 1st edge - 2nd node reported for the element
    edge_info((3*ielem)-1,1) = E2N(ielem,2);    %1st node of 2nd edge - 2nd node reported for the element
    edge_info((3*ielem)-1,2) = E2N(ielem,3);    %2nd node of 2nd edge - 3rd node reported for the element
    edge_info((3*ielem),1) = E2N(ielem,3);  %1st node of 3rd edge - 3rd node reported for the element
    edge_info((3*ielem),2) = E2N(ielem,1);  %2nd node of 3rd edge - 1st node reported for the element
    
    edge_info((3*ielem)-2, 3) = 3;
    edge_info((3*ielem)-1, 3) = 1;
    edge_info((3*ielem), 3) = 2;
    
    edge_info((3*ielem)-2,4) = ielem;    % Storing the element id for which the edge was determined
    edge_info((3*ielem)-1,4) = ielem;
    edge_info((3*ielem),4) = ielem;
end

% 4th column of edge infowill store the other bounding element, determined in the next steps
for i=1:n_edge_inc
    for j=1:nbfgrp
        check = any(l{j}==edge_info(i,1));
        if any(check > 0)    %checks if for any edge both nodes are listed in the list of nodes of a boundary face 
            edge_info(i,5) = -1*j;
        end
    end
end

ind = edge_info(:,5) >= 0;  % indices of edges which need to be filtered out 
edge_info(ind,:) =[];   %entries corresponding to those indices are removed

end