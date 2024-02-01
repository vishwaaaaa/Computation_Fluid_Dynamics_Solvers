function solnAdaptation(gri_in_file,fstatein,gamma,gri_out_file,fstateout)

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
  nnode=nnode_i(ibfgrp);
  for ibface = 1:nbface,
    Ash = fscanf(fid, '%d', nnode);
  end
  l{ibfgrp}=Ash;
  for i=1:height(Ash)
      node_b_check(Ash(i),ibfgrp)=1;
  end
end

% Read in elements and plot edges
figure(1); clf; hold on;
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


elem_state = zeros(nelem, 4);
elem_M = zeros(1,nelem);

fid2 = fopen(fstatein, 'r');
% Read in nodes

for i = 1:nelem,
  A2 = fscanf(fid2, '%lf', 4);
  elem_state(i,:) = A2(1:4)';
end

fclose(fid2);

for i=1:nelem
    p = (gamma-1)*(elem_state(i,4)-((elem_state(i,2)^2 + elem_state(i,3)^2)/(2*elem_state(i,1))));
    elem_M(i)=sqrt((elem_state(i,2)^2 + elem_state(i,3)^2)/(gamma*p*elem_state(i,1)));
end
%% Section 3 - Finding edges to be split

edge_error = zeros(1,n_edge_c);
edge_ref_c = 0;
for i=1:n_edge_c
    k = edge_info(i,4);
    if edge_info(i,4)>0
        edge_error(i) = abs(elem_M(edge_info(i,3))-elem_M(edge_info(i,4)))*sqrt(edge_info(i,5));
        edge_ref_c = edge_ref_c+1;
    else
        name = title{-1*k};
        name = string(name{1});
        if name =='Airfoil'|name =='airfoil'|name =='Flap'|name =='flap'|name =='Slat'|name =='slat'|name =='Main'|name =='Main'
            p = (gamma-1)*(elem_state(edge_info(i,3),4)-((elem_state(edge_info(i,3),2)^2 + elem_state(edge_info(i,3),3)^2)/(2*elem_state(edge_info(i,3),1))));
            c = sqrt(gamma*p/elem_state(edge_info(i,3),1));
            nx = (V(edge_info(i,2),2)-V(edge_info(i,1),2))/edge_info(i,5);
            ny = (V(edge_info(i,1),1)-V(edge_info(i,2),1))/edge_info(i,5);
            u_n = (elem_state(edge_info(i,3),2))*nx/elem_state(edge_info(i,3),1);
            v_n = (elem_state(edge_info(i,3),3))*nx/elem_state(edge_info(i,3),1);
            M = sqrt((u_n*u_n)+(v_n*v_n))/c;
            edge_error(i) =M*sqrt(edge_info(i,5));
            edge_ref_c = edge_ref_c+1;
        else 
            edge_error(i) = 0;
        end
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
elem_state_mod = repmat(elem_state,1);

for i=1:nelem
    if sum(E2IB_split_check(i,:))==3
        E3_1 = [E2N(i,1) edge_mp(E2IB(i,3)) edge_mp(E2IB(i,2))];
        E3_2 = [E2N(i,2) edge_mp(E2IB(i,1)) edge_mp(E2IB(i,3))];
        E3_3 = [E2N(i,3) edge_mp(E2IB(i,2)) edge_mp(E2IB(i,1))];
        E3_4 = [edge_mp(E2IB(i,1)) edge_mp(E2IB(i,2)) edge_mp(E2IB(i,3))];
        E2N_mod = [E2N_mod;E3_1;E3_2;E3_3;E3_4];
        elem_state_mod = [elem_state_mod;elem_state(i,:);elem_state(i,:);elem_state(i,:);elem_state(i,:)];
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
        elseif e1_id == 2
            E1_1 = [E2N(i,2) E2N(i,3) mp_req];
            E1_2 = [E2N(i,2) mp_req E2N(i,1)];
        elseif e1_id == 3
            E1_1 = [E2N(i,3) E2N(i,1) mp_req];
            E1_2 = [E2N(i,3) mp_req E2N(i,2)];
        end
        E2N_mod = [E2N_mod;E1_1;E1_2];
        elem_state_mod = [elem_state_mod;elem_state(i,:);elem_state(i,:)];
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
            E2_2 = [E2N(i,2) mp_req_1 mp_req_2];
            E2_3 = [E2N(i,3) mp_req_1 E2N(i,2)];
        elseif e2_id==2
            mp_req_1 = edge_mp(E2IB(i,3));
            mp_req_2 = edge_mp(E2IB(i,1));
            E2_1 = [E2N(i,2) mp_req_2 mp_req_1];
            E2_2 = [E2N(i,3) mp_req_1 mp_req_2];
            E2_3 = [E2N(i,1) mp_req_1 E2N(i,3)];
        elseif e2_id==3
            mp_req_1 = edge_mp(E2IB(i,1));
            mp_req_2 = edge_mp(E2IB(i,2));
            E2_1 = [E2N(i,3) mp_req_2 mp_req_1];
            E2_2 = [E2N(i,1) mp_req_1 mp_req_2];
            E2_3 = [E2N(i,2) mp_req_1 E2N(i,1)];
        end
        E2N_mod = [E2N_mod;E2_1;E2_2;E2_3];
        elem_state_mod = [elem_state_mod;elem_state(i,:);elem_state(i,:);elem_state(i,:)];
    end
end

r_e_id = [];
for i=1:nelem
    if sum(E2IB_split_check(i,:))~=0
        r_e_id = [r_e_id i];
    end
end

E2N_mod(r_e_id,:)=[];
elem_state_mod(r_e_id,:) = [];

%% Section 6 - Mentioning the new nodes in the boundary list

for i=1:n_edge_c
    if edge_info(i,4)<0 && edge_info(i,6)==1
        l{-1*edge_info(i,4)} = [l{-1*edge_info(i,4)};edge_mp(i)];
    end
end

nnode_b_new = zeros(nbfgrp,1);
for i=1:nbfgrp
    nnode_b_new(i) = length(l{i});
end

%% Section 7 - Spline Interpolation

airfoil_node_id_ref =[];
flap_node_id_ref =[];
slat_node_id_ref = [];
for ibfgrp = 1:nbfgrp
    name = title{ibfgrp};
    name = string(name{1});
    lt = l{ibfgrp};
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

%% Section 8 - output file writing

grioutf = fopen(gri_out_file,'w');
fprintf(grioutf,'%d %d %d\n',height(V_mod),height(E2N_mod),dim);
for i=1:length(V_mod)
    fprintf(grioutf,'%.10f %.10f\n',V_mod(i,1),V_mod(i,2));
end
fprintf(grioutf,'%d\n',nbfgrp);
for i=1:nbfgrp
    s = title{i};
    %lt = bnode_list_all{i};
    lt = l{i};
    fprintf(grioutf,'%d %d %s\n',nbface,nnode_b_new(i),s{1});
    for j=1:nnode_b_new(i)
        if j~=nnode_b_new(i) 
            fprintf(grioutf,'%d ',lt(j));
        else
            fprintf(grioutf,'%d\n',lt(j));
        end
        
    end
end
fprintf(grioutf,'%d %d %s',height(E2N_mod),p_in,sbasis{1});
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
