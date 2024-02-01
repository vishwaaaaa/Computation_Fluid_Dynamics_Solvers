function write_curve_gri(in_gri,node_loc,node_add, elem_rem, out_gri_name, Q)
%write_curve_gri creates the curved gri file from the original file
%INPUTS: in_gri is the input gri file
%        node_loc array of the new nodes added for the Q value
%        node_add array of the new nodes added to each element
%        out_gri_name is the name of what the new gri file should be
%        Q is the order of the curving

% reading in original data
fid = fopen(in_gri, 'r');
% Read in nodes
A = fscanf(fid,'%d', 3);
nnode    = A(1);
nelemtot = A(2);
dim      = A(3);
V = zeros(nnode, dim);
for inode = 1:nnode,
  A = fscanf(fid, '%lf', 2);
  V(inode,:) = A(1:2)';
end

A = fscanf(fid, '%d', 1);
nbfgrp = A(1);
for ibfgrp = 1:nbfgrp,
  fgets(fid);
  sline = fgets(fid);
  [nbface, nnode, title] = strread(sline, '%d %d %s');
  A = zeros(nbface, 2);
  for ibface = 1:nbface,
    Ash = fscanf(fid, '%d', nnode);
    A(ibface,:) = Ash;
  end
  l{ibfgrp}=A;
end

% Read in elements and plot edges
curtot = 0;
E2N = zeros(nelemtot, 3);
while (curtot ~= nelemtot),
  fgets(fid);
  sline = fgets(fid);
  [nelem, p_in, sbasis] = strread(sline, '%d %d %s');
  switch sbasis{1}
    case 'TriLagrange'
      nnode = (p_in+1)*(p_in+2)/2;
      nedge = 3;
      fvec = zeros(3,p_in+1);
      fvec(1,:) = [1:(p_in+1)];
      v = p_in+1; d = p_in;
      for k=1:(p_in+1), fvec(2,k) = v; v = v+d; d = d-1; end;
      v = 1; d = p_in+1;
      for k=1:(p_in+1), fvec(3,k) = v; v = v+d; d = d-1; end;
    otherwise
      error('element type not understood');
  end
  for elem = 1:nelem,
    A = fscanf(fid, '%d', nnode);  % HO nodes included
    E2N(elem,1)=A(1);
    E2N(elem,2)=A(2);
    E2N(elem,3)=A(3);
   
  end
  curtot = curtot + nelem;
end
fclose(fid);
nelem = nelemtot;

% Creating output file
grioutf = fopen(out_gri_name,'w');

fprintf(grioutf,'%d %d %d\n',length(V) + length(node_loc),length(E2N),2);
for i=1:length(V)
    fprintf(grioutf,'%.10f %.10f\n',V(i,1),V(i,2));
end
for i=1:length(node_loc)
    fprintf(grioutf,'%.10f %.10f\n',node_loc(i,1),node_loc(i,2));
end

title = {'farfield', 'slat', 'main', 'flap'};
str = string(title);
fprintf(grioutf,'%d\n',nbfgrp);
for i=1:nbfgrp
    s = str(i);
    lt = l{i};
    nbface = length(lt);
    fprintf(grioutf,'%d %d %s\n',nbface,2,s);
    for j=1:nbface
        fprintf(grioutf,'%d %d\n', lt(j,1), lt(j,2));
    end
end

fprintf(grioutf,'%d %d %s',height(E2N) - height(elem_rem),p_in,sbasis{1});
for i=1:length(E2N)
    val = ismember(i, elem_rem);
    if val == 1
        continue
    else
        fprintf(grioutf,'\n%d %d %d',E2N(i,1),E2N(i,2),E2N(i,3));
    end
end

fprintf(grioutf,'\n%d %d %s',height(node_add),Q,sbasis{1});
for i=1:length(node_add)
    if Q == 2
        fprintf(grioutf, '\n%d %d %d %d %d %d', node_add(i,1), node_add(i,2), node_add(i,3), node_add(i,4), node_add(i,5), node_add(i,6));
    end
    
    if Q == 3
        fprintf(grioutf, '\n%d %d %d %d %d %d %d %d %d %d', node_add(i,1), node_add(i,2), node_add(i,3), node_add(i,4), node_add(i,5), node_add(i,6), node_add(i,7), node_add(i,8), node_add(i,9), node_add(i,10));
    end
end

fclose(grioutf);
    

end

