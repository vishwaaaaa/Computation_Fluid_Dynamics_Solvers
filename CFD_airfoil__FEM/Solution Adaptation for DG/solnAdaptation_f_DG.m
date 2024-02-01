function solnAdaptation_f_DG(gri_in_file,soln_order,fstatein,edgeMin,gri_out_file,fstateout)

rank = floor((soln_order+1)*(soln_order+2)/2);

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
node_b_check = zeros(nnode,nbfgrp);
for ibfgrp = 1:nbfgrp,
  fgets(fid);
  sline = fgets(fid);
  [nbface, nnode_i(ibfgrp), title{ibfgrp}] = strread(sline, '%d %d %s');
  boundary_node = zeros(nbface,nnode_i(ibfgrp));
  boundary_node_after = zeros(2*nbface,nnode_i(ibfgrp));

  nnode=nnode_i(ibfgrp);
  for ibface = 1:nbface,
    Ash = fscanf(fid, '%d', nnode);
    for j=1:nnode_i(ibfgrp)
        boundary_node(ibface,j) = Ash(j);
        boundary_node_after(2*ibface,j) = Ash(j);
        node_b_check(Ash(j),ibfgrp)=1;
    end
  end
  l{ibfgrp}=boundary_node;
  l_mod{ibfgrp}=boundary_node_after;
end

% Read in elements and plot edges
figure(1); clf; hold on;
curtot = 0;
E2N_b1 = zeros(nelemtot, 3);
E2N_b2 = zeros(nelemtot, 6);
E2N_b3 = zeros(nelemtot, 10);

while (curtot ~= nelemtot),
  fgets(fid);
  sline = fgets(fid);
  [nelem, q_in, sbasis] = strread(sline, '%d %d %s');
  switch sbasis{1}
    case 'TriLagrange'
      nnode = (q_in+1)*(q_in+2)/2;
      nedge = 3;
      fvec = zeros(3,q_in+1);
      fvec(1,:) = [1:(q_in+1)];
      v = q_in+1; d = q_in;
      for k=1:(q_in+1), fvec(2,k) = v; v = v+d; d = d-1; end;
      v = 1; d = q_in+1;
      for k=1:(q_in+1), fvec(3,k) = v; v = v+d; d = d-1; end;
    otherwise
      error('element type not understood');
  end
  
  for elem = 1:nelem,
    A = fscanf(fid, '%d', nnode);  % HO nodes included
    if q_in==1
        for j=1:3
            E2N_b1(elem,j) = A(j);
        end
    end
    
    if q_in==2
        for j=1:6
            E2N_b2(elem,j) = A(j);
        end
    end

    if q_in==3
        for j=1:10
            E2N_b3(elem,j) = A(j);
        end
    end

    for edge=1:nedge,
      plotedge(V(A(fvec(edge,:)),:));
    end
    
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

E2N_b1((E2N_b1(:,1)==0),:) = [];
E2N_b2((E2N_b2(:,1)==0),:) = [];
E2N_b3((E2N_b3(:,1)==0),:) = [];

% Since we are only refining linear elements, E2N = E2N_b1

E2N = E2N_b1;

%% Section 1 - Getting all edge data

n_edge_inc = 3*nelem;
edge_info = zeros(n_edge_inc,6);
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

dup_edge = [];
IE = [];

ENE_inc = zeros(nelem,3);    %Stores neighboring elements
ENE_count = zeros(1,nelem);

E2IB_inc = zeros(nelem,3);  %Stores edge of each elements
E2IB_count = zeros(1,nelem);
E2IB_edge_ref = zeros(nelem,3);

for i=1:n_edge_inc
    edge_info(i,5) = sqrt((V(edge_info(i,1),1)-V(edge_info(i,2),1))^2 + (V(edge_info(i,1),2)-V(edge_info(i,2),2))^2);   %edge length calculation
        for j=i:n_edge_inc
            if edge_info(i,1)==edge_info(j,2) && edge_info(i,2)==edge_info(j,1) && j~=i
                edge_info(i,4) = edge_info(j,3);
                dup_edge = [dup_edge j];
                IE = [IE;edge_info(i,1) edge_info(i,2) edge_info(i,3) edge_info(i,4) edge_info(i,5) edge_info(i,6)];
                break
            end
        end
end
% How to do above loop in O(N) operations ?!?!

edge_info(dup_edge,:) =[];
n_edge_c = height(edge_info);

for i=1:n_edge_c
    if edge_info(i,4)==0
        for j=1:nbfgrp
            if node_b_check(edge_info(i,1),j)==1 && node_b_check(edge_info(i,2),j)==1
                edge_info(i,4) = -1*j;
            end
        end
    end
end


%Find edges and neighbouring elements of an element

for i=1:n_edge_c
    k1 = edge_info(i,3);
    k2 = edge_info(i,4);
    E2IB_count(k1)=E2IB_count(k1)+1;
    E2IB_inc(k1,E2IB_count(k1))=i;
    ENE_count(k1) = ENE_count(k1)+1;
    ENE_inc(k1,ENE_count(k1)) = k2;
    if k2>0
       E2IB_count(k2)=E2IB_count(k2)+1;
       E2IB_inc(k2,E2IB_count(k2))=i;
       ENE_count(k2) = ENE_count(k2)+1;
       ENE_inc(k2,ENE_count(k2)) = k1;
    end        
end

% currently the elements/edges are in haphazard order, need to order them in
% (a) counterclockwise direction, and (b) respecting order in ehcih nodes
% are represented in E2N

E2IB = repmat(E2IB_inc,1);
ENE = repmat(ENE_inc,1);

edge_order = zeros(nelem,3);

for i=1:nelem
    count = 1;
    for j=1:3
        q = E2N(i,j);
        for k=1:3
            w = E2IB_inc(i,k);
            if(q~=edge_info(w,1) && q~=edge_info(w,2))
                edge_order(i,k) = count;
                count = count+1;
                break
            end
        end
    end
    for j=1:3
        E2IB(i,edge_order(i,j))=E2IB_inc(i,j);
        ENE(i,edge_order(i,j))=ENE_inc(i,j);
    end
end

%splitting into interior edges and boundary edges

n_IE = height(IE);
n_BE = n_edge_c-n_IE;

BE = repmat(edge_info,1);
ie = [];
be = [];

for i=1:n_edge_c
    if edge_info(i,4)>0
        ie = [ie i];
    else
        be = [be i];
    end
end

BE(ie,:)=[];

%% Section 2 - Read state file and calculate M

elem_state = zeros(nelem, rank, 4);
elem_state_list = zeros(nelem*rank,4);
edge_M = zeros(1,n_edge_c);
edge_elemL = zeros(1,n_edge_c);
edge_elemR = zeros(1,n_edge_c);

fid2 = fopen(edgeMin, 'r');
% Read in nodes

for i = 1:n_edge_c
  interm = fscanf(fid2, '%lf',3);
  edge_elemL(i) = floor(interm(1)); edge_elemR(i) = floor(interm(2)); edge_M(i) = interm(3);
end

fclose(fid2);
fid3 = fopen(fstatein, 'r');
for i=1:nelem
    for j=1:rank
        interm = fscanf(fid3, '%lf', 4);
        for k=1:4
            elem_state(i,j,k)=interm(k);
            elem_state_list(((i-1)*rank)+j,k) =interm(k);
        end
    end
end
fclose(fid3);

%% Section 3 - Finding edges to be split

edge_error = zeros(1,n_edge_c);
edge_ref_c = 0;
for i=1:n_edge_c
    k1 = edge_info(:,3)==edge_elemL(i);
    k2 = edge_info(:,4)==edge_elemR(i);
    k_ind = find((k1 & k2)==1);
    edge_error(k_ind) =  edge_M(i);
    if edge_M(i)>0
        edge_ref_c = edge_ref_c+1;
    end
end

e_req = round(0.03*edge_ref_c);
if e_req==0
    [B,I] = maxk(edge_error,1);
else
    [B,I] = maxk(edge_error,e_req);
end

for i=1:length(I)
    edge_info(I(i),6)=1;
end

E2IB_split_check = zeros(nelem,3);
edge_info_itmd = repmat(edge_info,1);

for i=1:n_edge_c
    if(edge_info_itmd(i,6)==1)
        k1 = edge_info_itmd(i,3);
        k2 = edge_info_itmd(i,4);
        for j=1:3
            E2IB_split_check(k1,j)=1;
            if k2>0
                E2IB_split_check(k2,j)=1;
            end
        end
    end
end

for i=1:nelem
    for j=1:3
        if E2IB_split_check(i,j)==1
            edge_info(E2IB(i,j),6)=1;
        end
    end
end

for i=1:n_edge_c
    if edge_info(i,6)==1
        k1 = edge_info(i,3);
        k2 = edge_info(i,4);
        for j=1:3
            if E2IB(k1,j)==i
                E2IB_split_check(k1,j)=1;
            end
            if k2>0
                if E2IB(k2,j)==i
                    E2IB_split_check(k2,j)=1;
                end
            end
        end
    end
end

%% Section 4 - Creatig new nodes

V_mod = repmat(V,1);
nelem_new = 0;
nnode_new = 0;
nn_new=[];
nn_elem=[];
n_new=zeros(1,2);
edge_mp = zeros(n_edge_c,1);
count = 0;

for i=1:n_edge_c
    if edge_info(i,6)==1
        for j=1:2
            n_new(j) = (V(edge_info(i,1),j)+V(edge_info(i,2),j))/2;
        end
        V_mod = [V_mod;n_new];
        edge_mp(i) = height(V_mod);
    end
end

%% Section 5 - Creating new elements and replacing existing elements with them +defining same state for the new elements

E2N_mod = repmat(E2N,1);
%elem_state_mod = repmat(elem_state,1);
elem_state_mod = repmat(elem_state_list,1);

for i=1:nelem
    if sum(E2IB_split_check(i,:))==3
        E3_1 = [E2N(i,1) edge_mp(E2IB(i,3)) edge_mp(E2IB(i,2))];
        E3_3 = [edge_mp(E2IB(i,3)) E2N(i,2) edge_mp(E2IB(i,1))];
        E3_4 = [edge_mp(E2IB(i,2)) edge_mp(E2IB(i,1)) E2N(i,3)];
        E3_2 = [edge_mp(E2IB(i,2)) edge_mp(E2IB(i,3)) edge_mp(E2IB(i,1))];
        E2N_mod = [E2N_mod;E3_1;E3_2;E3_3;E3_4];
        s = [1,1,1];
        U_basis_3D = elem_state(i,:,:);
        U_basis = zeros(width(U_basis_3D),4);
        for j=1:width(U_basis_3D)
            for k=1:4
                U_basis(j,k) = U_basis_3D(1,j,k);
            end
        end
        [elem_state_1,elem_state_2,elem_state_3,elem_state_4] = Basis_nn(soln_order,s,U_basis);
        elem_state_mod = [elem_state_mod;elem_state_1;elem_state_2;elem_state_3;elem_state_4];
    elseif sum(E2IB_split_check(i,:))==1
        for j=1:3
            if E2IB_split_check(i,j)==1
                e1_id = j;
                mp_req = edge_mp(E2IB(i,j));
                break
            end
        end
        if e1_id ==1
            E1_1 = [E2N(i,1) E2N(i,2) mp_req];
            E1_2 = [E2N(i,1) mp_req E2N(i,3)];
            s = [1,0,0];
        elseif e1_id == 2
            E1_2 = [mp_req E2N(i,2) E2N(i,3)];
            E1_1 = [E2N(i,1) E2N(i,2) mp_req];
            s = [0,1,0];
        elseif e1_id == 3
            E1_1 = [E2N(i,1) mp_req E2N(i,3)];
            E1_2 = [mp_req E2N(i,2) E2N(i,3) ];
            s = [0,0,1];
        end
        E2N_mod = [E2N_mod;E1_1;E1_2];
        U_basis_3D = elem_state(i,:,:);
        U_basis = zeros(width(U_basis_3D),4);
        for j=1:width(U_basis_3D)
            for k=1:4
                U_basis(j,k) = U_basis_3D(1,j,k);
            end
        end
        [elem_state_1,elem_state_2,elem_state_3,elem_state_4] = Basis_nn(soln_order,s,U_basis);
        elem_state_mod = [elem_state_mod;elem_state_1;elem_state_2];
    elseif sum(E2IB_split_check(i,:))==2
        for j=1:3
            if E2IB_split_check(i,j)==0
                e2_id = j;
                break
            end
        end
        if e2_id==1
            mp_req_1 = edge_mp(E2IB(i,2));
            mp_req_2 = edge_mp(E2IB(i,3));
            E2_1 = [E2N(i,1) mp_req_2 mp_req_1];
            E2_2 = [mp_req_2 E2N(i,2) mp_req_1];
            E2_3 = [mp_req_1 E2N(i,2) E2N(i,3)];
            s = [0,1,1];
        elseif e2_id==2
            mp_req_1 = edge_mp(E2IB(i,3));
            mp_req_2 = edge_mp(E2IB(i,1));
            E2_3 = [mp_req_1 E2N(i,2) mp_req_2];
            E2_2 = [E2N(i,1) mp_req_1 mp_req_2];
            E2_1 = [E2N(i,1) mp_req_2 E2N(i,3)];
            s = [1,0,1];
        elseif e2_id==3
            mp_req_1 = edge_mp(E2IB(i,1));
            mp_req_2 = edge_mp(E2IB(i,2));
            E2_3 = [mp_req_2 mp_req_1 E2N(i,3)];
            E2_1 = [E2N(i,1) mp_req_1 mp_req_2];
            E2_2 = [E2N(i,1) E2N(i,2) mp_req_1];
            s = [1,1,0];
        end
        E2N_mod = [E2N_mod;E2_1;E2_2;E2_3];
        U_basis_3D = elem_state(i,:,:);
        U_basis = zeros(width(U_basis_3D),4);
        for j=1:width(U_basis_3D)
            for k=1:4
                U_basis(j,k) = U_basis_3D(1,j,k);
            end
        end
        [elem_state_1,elem_state_2,elem_state_3,elem_state_4] = Basis_nn(soln_order,s,U_basis);
        elem_state_mod = [elem_state_mod;elem_state_1;elem_state_2;elem_state_3];
    end
end


r_e_id = [];
for i=1:nelem
    if sum(E2IB_split_check(i,:))~=0
        r_e_id = [r_e_id i];
    end
end

r_b_id = [];
for i=1:length(r_e_id)
    id1 = (rank*(r_e_id(i)-1))+1;
    id2 = rank*r_e_id(i);
    for j=id1:id2
        r_b_id = [r_b_id j];
    end
end

E2N_mod(r_e_id,:)=[];
elem_state_mod(r_b_id,:) = [];

%% Section 6 - Mentioning the new nodes in the boundary list

for i=1:n_edge_c
    if edge_info(i,4)<0 && edge_info(i,6)==1
        k1 = find(l_mod{-1*edge_info(i,4)}(:,1)==edge_info(i,1));
        if l_mod{-1*edge_info(i,4)}(k1,2)~=edge_info(i,2)
            k1 = find(l_mod{-1*edge_info(i,4)}(:,2)==edge_info(i,1));
        end
        l_mod{-1*edge_info(i,4)}(k1-1,1)=l_mod{-1*edge_info(i,4)}(k1,1);
        l_mod{-1*edge_info(i,4)}(k1-1,2)= edge_mp(i);
        l_mod{-1*edge_info(i,4)}(k1,1)= edge_mp(i);
    end
end

for i=1:nbfgrp
    Aqua = l_mod{i};
    Aqua((Aqua(:,1)==0),:) = [];
    l_mod{i} = Aqua;
end

nbface_new = zeros(nbfgrp,1);
nnode_b_new = zeros(nbfgrp,1);
for i=1:nbfgrp
    nnode_b_new(i) = width(l_mod{i});
    nbface_new(i) = height(l_mod{i});
end


%% Section 7 - Spline Interpolation

airfoil_node_id_ref =[];
flap_node_id_ref =[];
slat_node_id_ref = [];

for ibfgrp = 1:nbfgrp
    name = title{ibfgrp};
    name = string(name{1});
    %lt = l{ibfgrp};
    lt = transpose(l_mod{ibfgrp}(:,1));
    llt = length(lt);
    if name =='Airfoil'|name =='airfoil'
        for j=1:llt
            n=lt(j);
            V_mod(n,:) = nspline([V_mod(n,:)],2);
        end
    elseif name =='Flap'|name =='flap'
        for j=1:llt
            n=lt(j);
            V_mod(n,:) = nspline([V_mod(n,:)],1);
        end
    elseif name =='Slat'|name =='slat'
        for j=1:llt
            n=lt(j);
            V_mod(n,:) = nspline([V_mod(n,:)],3);
        end
    end
end
%% Section 8 - output file writing

grioutf = fopen(gri_out_file,'w');
fprintf(grioutf,'%d %d %d\n',height(V_mod),height(E2N_mod),dim);
for i=1:length(V_mod)
    fprintf(grioutf,'%.10f %.10f\n',V_mod(i,1),V_mod(i,2));
end
fprintf(grioutf,'%d\n',nbfgrp);
for i=1:nbfgrp
    s = title{i};
    lt = l_mod{i};
    fprintf(grioutf,'%d %d %s\n',nbface_new(i),nnode_b_new(i),s{1});
    for j=1:nbface_new(i)
        fprintf(grioutf,'%d',lt(j,1));
        for a=2:nnode_b_new(i)
            fprintf(grioutf,' %d',lt(j,a));
        end
        fprintf(grioutf,'\n');
    end
    
end
fprintf(grioutf,'%d %d %s',height(E2N_mod),q_in,sbasis{1});
for i=1:length(E2N_mod)
    fprintf(grioutf,'\n%d %d %d',E2N_mod(i,1),E2N_mod(i,2),E2N_mod(i,3));
end
fclose(grioutf);


state_out = fopen(fstateout,'w');
for i=1:length(elem_state_mod)
    fprintf(state_out,'%d %d %d %d',elem_state_mod(i,1),elem_state_mod(i,2),elem_state_mod(i,3),elem_state_mod(i,4));
    if i~=length(elem_state_mod)
        fprintf(state_out,'\n');
    end
end

fclose(state_out);

end

function plotedge(V)
x = V(:,1);
y = V(:,2);
plot(x,y, 'k-', 'linewidth', 1);
end
