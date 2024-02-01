function ref_uniform(gri_in_file,gri1_out_file)

fid = fopen(gri_in_file, 'r');
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
% Read through boundary info
A = fscanf(fid, '%d', 1);
nbfgrp = A(1);
nnode_i = zeros(nbfgrp,1);
for ibfgrp = 1:nbfgrp,
  fgets(fid);
  sline = fgets(fid);
  [nbface, nnode_i(ibfgrp), title{ibfgrp}] = strread(sline, '%d %d %s');
  nnode=nnode_i(ibfgrp);
  for ibface = 1:nbface,
    A = fscanf(fid, '%d', nnode);
  end
  l{ibfgrp}=A;
end

% Read in elements and plot edges
%figure(1); clf; hold on;
curtot = 0;
E2N = zeros(nelemtot, 3);
while (curtot ~= nelemtot),
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
      for k=1:(p+1), fvec(2,k) = v; v = v+d; d = d-1; end;
      v = 1; d = p+1;
      for k=1:(p+1), fvec(3,k) = v; v = v+d; d = d-1; end;
    otherwise
      error('element type not understood');
  end
  
  for elem = 1:nelem,
    A = fscanf(fid, '%d', nnode);  % HO nodes included
    E2N(elem,1)=A(1);
    E2N(elem,2)=A(2);
    E2N(elem,3)=A(3);
    %{
    for edge=1:nedge,
      plotedge(V(A(fvec(edge,:)),:));
    end
    %} 
  end
  curtot = curtot + nelem;
end
fclose(fid);
nelem = nelemtot;

axis equal
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
box off; axis off
set(gca, 'LineWidth', 0.5);

n_new=zeros(nelem*3,2);

count = 1;
for i=1:nelem
    for j=1:2
        n_new(count,j)= (V(E2N(i,1),j)+V(E2N(i,2),j))/2;
        n_new(count+1,j)= (V(E2N(i,2),j)+V(E2N(i,3),j))/2;
        n_new(count+2,j)= (V(E2N(i,3),j)+V(E2N(i,1),j))/2;
    end
    count = count+3;
end

n_new_c = unique(n_new(:,1:2),'rows','stable');
n_tot = [V; n_new_c];

nnode_i_new = (2*nnode_i)-1;
for ibfgrp = 1:nbfgrp
    bnode_list = zeros(nnode_i_new(ibfgrp),1);
    dayum = l{ibfgrp};
    for i=1:nnode_i(ibfgrp)
        bnode_list((2*i)-1)=dayum(i);
    end
    for i=1:nnode_i(ibfgrp)-1
        xr = (V(bnode_list((2*i)-1),1)+V(bnode_list((2*i)+1),1))/2;
        yr = (V(bnode_list((2*i)-1),2)+V(bnode_list((2*i)+1),2))/2;
        if(isempty(find(n_tot(:,1)==xr & n_tot(:,2)==yr))) %Changed this section to match ref_local and it works now
            bnode_list(2*i) =-10000;
        else
            bnode_list(2*i)=find(n_tot(:,1)==xr & n_tot(:,2)==yr);
        end
    end
    ind3 = bnode_list == -10000;
    bnode_list(ind3,:) =[];
    nnode_i_new(ibfgrp) = length(bnode_list);
    bnode_list_all{ibfgrp}=bnode_list;
end

nelemtot_new = nelemtot*4;
E2N_new = zeros(nelemtot_new,3);
count = 1;
for i=1:nelemtot
    N1 = E2N(i,1);
    N2 = E2N(i,2);
    N3 = E2N(i,3);
    n12 = find(n_tot(:,1)==(V(N1,1)+V(N2,1))/2 & n_tot(:,2)==(V(N1,2)+V(N2,2))/2);
    n23 = find(n_tot(:,1)==(V(N2,1)+V(N3,1))/2 & n_tot(:,2)==(V(N2,2)+V(N3,2))/2);
    n31 = find(n_tot(:,1)==(V(N3,1)+V(N1,1))/2 & n_tot(:,2)==(V(N3,2)+V(N1,2))/2);
    E2N_new(count,:)=[N1 n12 n31];
    E2N_new(count+1,:)=[N2 n23 n12];
    E2N_new(count+2,:)=[N3 n31 n23];
    E2N_new(count+3,:)=[n12 n23 n31];
    count = count+4;
end

for ibfgrp = 1:nbfgrp
    name = title{ibfgrp};
    name = string(name{1});
    lt = bnode_list_all{ibfgrp};
    llt = length(lt);
    if name =='Airfoil'|name =='airfoil'
        %airfoil_node = zeros(llt,2);
        for j=1:llt
            n=lt(j);
            %airfoil_node(j,:) = [V_mod(n,1) V_mod(n,2)];
            n_tot(n,:) = nspline([n_tot(n,:)],2);
        end
        %airfoil_node_ref = nspline(airfoil_node,2);
        %airfoil_node_id_ref = [lt airfoil_node_ref];
    elseif name =='Flap'|name =='flap'
%         flap_node = zeros(llt,2);
        for j=1:llt
            n=lt(j);
            %flap_node(j,:) = [V_mod(n,1) V_mod(n,2)];
            n_tot(n,:) = nspline([n_tot(n,:)],1);
        end
%         flap_node_ref = nspline(flap_node,1);
%         flap_node_id_ref = [lt flap_node_ref];
    elseif name =='Slat'|name =='slat'
%         slat_node = zeros(llt,2);
        for j=1:llt
            n=lt(j);
%             slat_node(j,:) = [V_mod(n,1) V_mod(n,2)];
            n_tot(n,:) = nspline([n_tot(n,:)],3);
        end
%    `    slat_node_ref = nspline(slat_node,3);
%         slat_node_id_ref = [lt slat_node_ref];
    end
end

% modified_nodes = [airfoil_node_id_ref;flap_node_id_ref;slat_node_id_ref];
% 
% if isempty(modified_nodes)
% else
%     for i=1:length(modified_nodes)
%         n = modified_nodes(i,1);
%         V_mod(n,1) = modified_nodes(i,2);
%         V_mod(n,1) = modified_nodes(i,3);
%     end
% end

gri_op = fopen(gri1_out_file,'w');
fprintf(gri_op,'%d %d %d\n',length(n_tot),length(E2N_new),dim);
for i=1:length(n_tot)
    %fprintf(gri_op,'%lf %lf\n',n_tot(i,1),n_tot(i,2));
    fprintf(gri_op,'%.10f %.10f\n',n_tot(i,1),n_tot(i,2));
end
fprintf(gri_op,'%d\n',nbfgrp);
for i=1:nbfgrp
    s = title{i};
    lt = bnode_list_all{i};
    fprintf(gri_op,'%d %d %s\n',nbface,nnode_i_new(i),s{1});
    for j=1:nnode_i_new(i)
        if j~=nnode_i_new(i) 
            fprintf(gri_op,'%d ',lt(j));
        else
            fprintf(gri_op,'%d\n',lt(j));
        end
        
    end
end
fprintf(gri_op,'%d %d %s',nelemtot_new,p,sbasis{1});
for i=1:nelemtot_new
    fprintf(gri_op,'\n%d %d %d',E2N_new(i,1),E2N_new(i,2),E2N_new(i,3));
end
fclose(gri_op);

end

function plotedge(V)
x = V(:,1);
y = V(:,2);
plot(x,y, 'k-', 'linewidth', 1);
end

