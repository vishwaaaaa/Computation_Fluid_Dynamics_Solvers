function ref_local(gri_in_file,x_lr,y_lr,r_lr,omega, gri1_out_file)

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
figure(1); clf; hold on;
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
%% I start from here 

n_edge_inc = 3*nelem;
edge_info = zeros(n_edge_inc,4);
for ielem=1:nelem
    edge_info((3*ielem)-2,1) = E2N(ielem,1);
    edge_info((3*ielem)-2,2) = E2N(ielem,2);
    edge_info((3*ielem)-1,1) = E2N(ielem,2);
    edge_info((3*ielem)-1,2) = E2N(ielem,3);
    edge_info((3*ielem),1) = E2N(ielem,3);
    edge_info((3*ielem),2) = E2N(ielem,1);
    edge_info((3*ielem)-2,3) = ielem;
    edge_info((3*ielem)-1,3) = ielem;
    edge_info((3*ielem),3) = ielem;
end

for i=1:n_edge_inc
    for j=1:nbfgrp
        if any(l{j}==edge_info(i,1)) && any(l{j}==edge_info(i,2))
            edge_info(i,4) = -1*j;
        end
    end
end

for i=1:n_edge_inc
    if edge_info(i,4)==0
        for j=i:n_edge_inc
            if edge_info(i,1)==edge_info(j,2) && edge_info(i,2)==edge_info(j,1) && j~=i
                edge_info(i,4) = edge_info(j,3);
                %edge_info(j,:) = [];
                edge_info(j,4) = -10000;
                break
            end
        end
    end
end

ind = edge_info(:,4) == -10000;
edge_info(ind,:) =[];
n_edge_c = length(edge_info);

ref_ch_elem = zeros(nelem,1);
for i=1:nelem
    x = [V(E2N(i,1),1) V(E2N(i,2),1) V(E2N(i,3),1)];
    y = [V(E2N(i,1),2) V(E2N(i,2),2) V(E2N(i,3),2)];
    polyin = polyshape({x},{y});
    [cen_x,cen_y] = centroid(polyin);
    d_cen_lrp = sqrt(((cen_x-x_lr)^2)+((cen_y-y_lr)^2));
    if d_cen_lrp>r_lr
       ref_ch_elem(i)=0;
    else 
       ref_ch_elem(i) =3;
    end
end

ref_ch_edge = zeros(n_edge_c,1);
for i=1:n_edge_c
    if edge_info(i,4)<0 && ref_ch_elem(edge_info(i,3))==3
        ref_ch_edge(i)=1;
    elseif edge_info(i,4) >0 && ref_ch_elem(edge_info(i,3))==3 && ref_ch_elem(edge_info(i,4))<3
        ref_ch_edge(i)=1;
        ref_ch_elem(edge_info(i,4)) = ref_ch_elem(edge_info(i,4))+1;
    elseif edge_info(i,4) >0 && ref_ch_elem(edge_info(i,3))<3 && ref_ch_elem(edge_info(i,4))==3
        ref_ch_edge(i)=1;
        ref_ch_elem(edge_info(i,3)) = ref_ch_elem(edge_info(i,3))+1;
    elseif edge_info(i,4) >0 && ref_ch_elem(edge_info(i,3))==3 && ref_ch_elem(edge_info(i,4))==3
        ref_ch_edge(i)=1;
    end   
end

E2N_mod = repmat(E2N,1);
V_mod = repmat(V,1);
E2N_new = [];

nelem_new = 0;
nnode_new = 0;
nn_new=[];
nn_elem=[];
for ielem=1:nelem
    if ref_ch_elem(ielem)==3
        n_new=zeros(3,2);
        count=1;
        for j=1:2
            n_new(count,j)= (V(E2N(ielem,1),j)+V(E2N(ielem,2),j))/2;
            n_new(count+1,j)= (V(E2N(ielem,2),j)+V(E2N(ielem,3),j))/2;
            n_new(count+2,j)= (V(E2N(ielem,3),j)+V(E2N(ielem,1),j))/2;
        end
        nn_new=[nn_new;n_new];
        nnode_new = nnode_new+3;
    end
end


V_mod = [V_mod;nn_new];
V_mod = unique(V_mod(:,1:2),'rows','stable');
edge_inf_all = [edge_info ref_ch_edge];



for i=1:nelem
    if ref_ch_elem(i)==3
        N1 = E2N(i,1);
        N2 = E2N(i,2);
        N3 = E2N(i,3);
        n12 = find(V_mod(:,1)==(V(N1,1)+V(N2,1))/2 & V_mod(:,2)==(V(N1,2)+V(N2,2))/2);
        n23 = find(V_mod(:,1)==(V(N2,1)+V(N3,1))/2 & V_mod(:,2)==(V(N2,2)+V(N3,2))/2);
        n31 = find(V_mod(:,1)==(V(N3,1)+V(N1,1))/2 & V_mod(:,2)==(V(N3,2)+V(N1,2))/2);
        E2N_new = [E2N_new;[N1 n12 n31];[N2 n23 n12];[N3 n31 n23];[n12 n23 n31]];
        nelem_new = nelem_new + 4;
    elseif ref_ch_elem(i)==1
        j1 = find(edge_inf_all(:,3)==i & edge_inf_all(:,5)==1);
        j2 = find(edge_inf_all(:,4)==i & edge_inf_all(:,5)==1);
        if(isempty(j1))
            e = j2;
            N1 = edge_inf_all(e,1);
            N2 = edge_inf_all(e,2);
            n12 = find(V_mod(:,1)==(V(N1,1)+V(N2,1))/2 & V_mod(:,2)==(V(N1,2)+V(N2,2))/2);
            bo = find(E2N(i,:)~=N1 & E2N(i,:)~=N2);
            N3 = E2N(i,bo);
            E2N_new = [E2N_new;[N2 n12 N3];[N3 n12 N1]];
        elseif(isempty(j2))
            e = j1;
            N1 = edge_inf_all(e,1);
            N2 = edge_inf_all(e,2);
            n12 = find(V_mod(:,1)==(V(N1,1)+V(N2,1))/2 & V_mod(:,2)==(V(N1,2)+V(N2,2))/2);
            bo = find(E2N(i,:)~=N1 & E2N(i,:)~=N2);
            N3 = E2N(i,bo);
            E2N_new = [E2N_new;[N1 n12 N3];[N3 n12 N2]];
        end
        nelem_new = nelem_new + 2;
    elseif ref_ch_elem(i)==2
        e = find((edge_inf_all(:,3)==i|edge_inf_all(:,4)==i) & edge_inf_all(:,5)==1);
        N_1 = edge_inf_all(e(1),1);
        N_2 = edge_inf_all(e(1),2);
        N_3 = edge_inf_all(e(2),1);
        N_4 = edge_inf_all(e(2),2);
        NL = [N_1 N_2 N_3 N_4];
        %Nl = E2N(i,:);
        %Ncount = histcounts(NL, Nl);
        Nl = unique(NL);
        Ncount = histcounts(NL, Nl);
        NR = find(Ncount==2);
        if(NR==1)
            nM2 = find(V_mod(:,1)==(V(Nl(3),1)+V(Nl(1),1))/2 & V_mod(:,2)==(V(Nl(3),2)+V(Nl(1),2))/2);
            nM1 = find(V_mod(:,1)==(V(Nl(1),1)+V(Nl(2),1))/2 & V_mod(:,2)==(V(Nl(1),2)+V(Nl(2),2))/2);
            E2N_new = [E2N_new;[Nl(1) nM1 nM2];[nM2 nM1 Nl(3)];[Nl(3) nM1 Nl(2)]];
        elseif(NR==2)
            nM2 = find(V_mod(:,1)==(V(Nl(1),1)+V(Nl(2),1))/2 & V_mod(:,2)==(V(Nl(1),2)+V(Nl(2),2))/2);
            nM1 = find(V_mod(:,1)==(V(Nl(3),1)+V(Nl(2),1))/2 & V_mod(:,2)==(V(Nl(3),2)+V(Nl(2),2))/2);
            E2N_new = [E2N_new;[Nl(2) nM1 nM2];[nM2 nM1 Nl(1)];[Nl(1) nM1 Nl(3)]];
        elseif(NR==3)
            nM1 = find(V_mod(:,1)==(V(Nl(3),1)+V(Nl(1),1))/2 & V_mod(:,2)==(V(Nl(3),2)+V(Nl(1),2))/2);
            nM2 = find(V_mod(:,1)==(V(Nl(3),1)+V(Nl(2),1))/2 & V_mod(:,2)==(V(Nl(3),2)+V(Nl(2),2))/2);
            E2N_new = [E2N_new;[Nl(3) nM1 nM2];[nM2 nM1 Nl(2)];[Nl(2) nM1 Nl(1)]];
        end
        nelem_new = nelem_new + 3;
        end
    end

ind2 = ref_ch_elem(:) > 0;
E2N_mod(ind2,:) =[];
E2N_mod = [E2N_mod;E2N_new];

nn_new = unique(nn_new(:,1:2),'rows','stable');
adrak = transpose(linspace(length(V)+1,length(V_mod),length(V_mod)-length(V)));
nn_new = [adrak nn_new];
old_presmoothing = repmat(edge_inf_all,1);

c1 = find(old_presmoothing(:,4)<0);
old_presmoothing(c1,:)=[];
c2 = find(old_presmoothing(:,5)==0);
old_presmoothing(c2,:)=[];
nn_old_presmoothing = [old_presmoothing(:,1);old_presmoothing(:,2)];
nn_old_presmoothing = unique(nn_old_presmoothing,'rows');

nn_new_presmooth = repmat(nn_new,1);

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
        if(isempty(find(V_mod(:,1)==xr & V_mod(:,2)==yr)))
            bnode_list(2*i) =-10000;
        else
            bnode_list(2*i)=find(V_mod(:,1)==xr & V_mod(:,2)==yr);
        end
    end
    ind3 = bnode_list == -10000;
    bnode_list(ind3,:) =[];
    nnode_i_new(ibfgrp) = length(bnode_list);
    bnode_list_all{ibfgrp}=bnode_list;
    for i=1:length(bnode_list)
        ind4 = find(nn_new_presmooth(:,1)==bnode_list(i));
        ind5 = find(nn_old_presmoothing(:,1)==bnode_list(i));
        nn_new_presmooth(ind4,:)=[];
        nn_old_presmoothing(ind5,:)=[];
    end
end

nn_old_presmooth=[];
for i=1:length(nn_old_presmoothing)
    nn_old_presmooth(i,:) = [nn_old_presmoothing(i,:) V(nn_old_presmoothing(i,:),1) V(nn_old_presmoothing(i,:),2)];
end
nn_presmooth = [nn_old_presmooth;nn_new_presmooth];
nn_smooth = repmat(nn_presmooth,1);

for i=1:length(nn_presmooth(:,1))
    smoothing = repmat(nn_presmooth,1);
    smoothing(i,:)=[];
    Ns = length(smoothing(:,1)); 
    X_sum = sum(smoothing(:,2));
    Y_sum = sum(smoothing(:,3));
    if(Ns ~=0)
        nn_smooth(i,2)=((1-omega)*nn_presmooth(i,2))+(omega*X_sum/Ns);
        nn_smooth(i,3)=((1-omega)*nn_presmooth(i,3))+(omega*Y_sum/Ns);
    end
end

for i=1:length(nn_smooth(:,1))
    j=nn_smooth(i);
    V_mod(j,1)=nn_smooth(i,2);
    V_mod(j,2)=nn_smooth(i,3);
    %for j=1:length(V_mod)
    %    if(nn_smooth(i)==j)
    %        
    %        break;
    %    end
    %end
end
airfoil_node_id_ref =[];
flap_node_id_ref =[];
slat_node_id_ref = [];
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
            V_mod(n,:) = nspline([V_mod(n,:)],2);
        end
        %airfoil_node_ref = nspline(airfoil_node,2);
        %airfoil_node_id_ref = [lt airfoil_node_ref];
    elseif name =='Flap'|name =='flap'
%         flap_node = zeros(llt,2);
        for j=1:llt
            n=lt(j);
            %flap_node(j,:) = [V_mod(n,1) V_mod(n,2)];
            V_mod(n,:) = nspline([V_mod(n,:)],1);
        end
%         flap_node_ref = nspline(flap_node,1);
%         flap_node_id_ref = [lt flap_node_ref];
    elseif name =='Slat'|name =='slat'
%         slat_node = zeros(llt,2);
        for j=1:llt
            n=lt(j);
%             slat_node(j,:) = [V_mod(n,1) V_mod(n,2)];
            V_mod(n,:) = nspline([V_mod(n,:)],3);
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

grioutf = fopen(gri1_out_file,'w');
fprintf(grioutf,'%d %d %d\n',length(V_mod),length(E2N_mod),dim);
for i=1:length(V_mod)
    fprintf(grioutf,'%.10f %.10f\n',V_mod(i,1),V_mod(i,2));
end
fprintf(grioutf,'%d\n',nbfgrp);
for i=1:nbfgrp
    s = title{i};
    lt = bnode_list_all{i};
    fprintf(grioutf,'%d %d %s\n',nbface,nnode_i_new(i),s{1});
    for j=1:nnode_i_new(i)
        if j~=nnode_i_new(i) 
            fprintf(grioutf,'%d ',lt(j));
        else
            fprintf(grioutf,'%d\n',lt(j));
        end
        
    end
end
fprintf(grioutf,'%d %d %s',length(E2N_mod),p,sbasis{1});
for i=1:length(E2N_mod)
    fprintf(grioutf,'\n%d %d %d',E2N_mod(i,1),E2N_mod(i,2),E2N_mod(i,3));
end
fclose(grioutf);
end

function plotedge(V)
x = V(:,1);
y = V(:,2);
plot(x,y, 'k-', 'linewidth', 1);
end