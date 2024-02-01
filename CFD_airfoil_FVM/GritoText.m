%AERO623 
%Make .txt files of essential info (nodes, elems, connect) from .gri
clear all; close all; clc;

fid = fopen("C:\Users\benba\Documents\MATLAB\Adaptive Refinement Function\Adaptive Refinement Function\c0_ED_02_It5.gri", 'r');
% Read in nodes
A = fscanf(fid,'%d', 3);
nnode    = A(1);
nelemtot = A(2);
dim      = A(3);
N = zeros(nnode, dim);
for inode = 1:nnode
  A = fscanf(fid, '%lf', 2);
  N(inode,:) = A(1:2)';
end
% Read through boundary info
A = fscanf(fid, '%d', 1);
nbfgrp = A(1);
NB = 0; NA = 0;
for ibfgrp = 1:nbfgrp
  fgets(fid);
  sline = fgets(fid);
  [nbface, nnode, title] = strread(sline, '%d %d %s');
  nnode;
  %for ibface = 1:nnode
    A = fscanf(fid, '%d', nnode);
    if ibfgrp < 2
        NB = [NB; A];
    else
        NA = [NA; A];
    end
  %end
end
NB = NB(2:end); NA = NA(2:end);
% Read through elements
fgets(fid);
line = fgets(fid);
line = textscan(line, '%d %d %s');
nelem = line{1};
E = zeros(nelem,3);
for i=1:nelem
    line = fgets(fid);
    line = textscan(line,'%d %d %d');
    E(i,:) = cell2mat(line);
end
[IE,BE] = edgehash(E);
C = zeros(nelem,3);
for i=1:height(IE)
    e1 = IE(i,3); e2 = IE(i,4); 
    [~,n1,~] = intersect(E(e1,:),IE(i,1:2));
    ii = ones(1,3); ii(n1) = 0; ii = find(ii~=0);
    C(e1,ii) = e2;
    [~,n2,~] = intersect(E(e2,:),IE(i,1:2));
    ii = ones(1,3); ii(n2) = 0; ii = find(ii~=0);
    C(e2,ii) = e1;
end
for i=1:height(BE)
    if any(NA == BE(i,1))
        e1 = BE(i,3);
        [~,n1,~] = intersect(E(e1,:),BE(i,1:2));
        ii = ones(1,3); ii(n1) = 0; ii = find(ii~=0);
        C(e1,ii) = -2;
    else
        e1 = BE(i,3);
        [~,n1,~] = intersect(E(e1,:),BE(i,1:2));
        ii = ones(1,3); ii(n1) = 0; ii = find(ii~=0);
        C(e1,ii) = -1;
    end
end

fidN = fopen("C:\Users\benba\Documents\C++\AERO623Proj2\SecondOrderAdapt5Nodes.txt",'w');
for i=1:height(N)
    fprintf(fidN,'%1.4f %1.4f \r\n',N(i,:));
end
fidE = fopen("C:\Users\benba\Documents\C++\AERO623Proj2\SecondOrderAdapt5Elems.txt",'w');
for i=1:height(E)
    fprintf(fidE,'%i %i %i \r\n',E(i,:));
end
fidC = fopen("C:\Users\benba\Documents\C++\AERO623Proj2\SecondOrderAdapt5Connect.txt",'w');
for i=1:height(C)
    fprintf(fidC,'%i %i %i \r\n',C(i,:));
end


function [IE,BE] = edgehash(E2N)
% This function identifies interior and boundary edges, and their
% connectivities, in a triangular mesh given an element-to-node array.
%
% INPUT : E2N = [nelem x 3] array mapping triangles to nodes
% OUTPUT: IE = [niedge x 4] array giving (n1, n2, elem1, elem2)
%              information for each interior edge
%         BE = [nbedge x 3] array giving (n1, n2, elem)
%              information for each boundary edge
    
    nelem = size(E2N,1);             % number of elements
    nnode = max(max(E2N));           % number of nodes
    H = sparse(nnode,nnode);        % Create a hash list to identify edges
    IE = zeros(ceil(nelem*3/2), 4); % (over) allocate interior edge array
    niedge = 0;                     % number of interior edges (running total)

    % Loop over elements and identify all edges
    for elem = 1:nelem
        nv = E2N(elem,1:3);
        for edge = 1:3
            n1 = nv(mod(edge  ,3)+1);
            n2 = nv(mod(edge+1,3)+1);
            if (H(n1,n2) == 0) % edge hit for the first time
                % could be a boundary or interior; assume boundary
                H(n1,n2) = elem;  H(n2,n1) = elem;
            else % this is an interior edge, hit for the second time
                oldelem = H(n1,n2);
                if (oldelem < 0), error 'Mesh input error'; end
                niedge = niedge+1;
                IE(niedge,:) = [n1,n2, oldelem, elem];
                H(n1,n2) = -1;  H(n2,n1) = -1;
            end
        end
    end

    IE = IE(1:niedge,:);  % clip IE

    % find boundary edges
    [I,J] = find(triu(H)>0);
    BE = [I, J, zeros(size(I))];
    for b = 1:size(I,1), BE(b, 3) = H(I(b),J(b));
    end
end
