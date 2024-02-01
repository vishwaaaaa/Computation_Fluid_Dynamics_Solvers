function makegri(N,E)
    E = triccw(E,N);
    nNode = height(N);
    x = N(:,1);
    y = N(:,2);
    nElemTot = height(E);
    Dim = 2;
    nBGroup = 7;
    nBFace = 1;
    Title = ["Bottom" "Right" "Top" "Left" "Slat" "Main" "Flap"];
    NB = zeros(7,200,1);
    xmin = min(N(:,1));
    xmax = max(N(:,1));
    ymax = max(N(:,2));
    ymin = min(N(:,2));
    x1L = -0.5;
    x1R = -0.01;
    x2L = x1R;
    x2R = 1.01;
    x3L = x2R;
    x3R = 1.5;
    [~,SB]=edgehash(E);
    t = find((N(SB(:,1),1)<xmin+1e-3)&(N(SB(:,2),1)<xmin+1e-3));
    t = unique(SB(t,1:2));
    nf(4) = length(t);
    NB(4,1:length(t)) = t;
    t = find((N(SB(:,1),1)>xmax-1e-3)&(N(SB(:,2),1)>xmax-1e-3));
    t = unique(SB(t,1:2));
    nf(2) = length(t);
    NB(2,1:length(t),:) = t;
    t = find((N(SB(:,1),2)<ymin+1e-3)&(N(SB(:,2),2)<ymin+1e-3));
    t = unique(SB(t,1:2));
    nf(1) = length(t);
    NB(1,1:length(t),:) = t;
    t = find((N(SB(:,1),2)>ymax-1e-3)&(N(SB(:,2),2)>ymax-1e-3));
    t = unique(SB(t,1:2));
    nf(3) = length(t);
    NB(3,1:length(t),:) = t;
    t1 = (N(SB(:,1),1)>x1L)&(N(SB(:,2),1)>x1L);
    t2 = (N(SB(:,1),1)<x1R)&(N(SB(:,2),1)<x1R);
    t = find(t1&t2);
    t = unique(SB(t,1:2));
    nf(5) = length(t);
    NB(5,1:length(t),:) = t;
    t1 = (N(SB(:,1),1)>x2L)&(N(SB(:,2),1)>x2L);
    t2 = (N(SB(:,1),1)<x2R)&(N(SB(:,2),1)<x2R);
    t = find(t1&t2);
    t = unique(SB(t,1:2));
    nf(6) = length(t);
    NB(6,1:length(t),:) = t;
    t1 = (N(SB(:,1),1)>x3L)&(N(SB(:,2),1)>x3L);
    t2 = (N(SB(:,1),1)<x3R)&(N(SB(:,2),1)<x3R);
    t = find(t1&t2);
    t = unique(SB(t,1:2));
    nf(7) = length(t);
    NB(7,1:length(t),:) = t;
    nElemGroup = 1;
    nElem = nElemTot;
    Order = 1;
    Basis = "TriLagrange";
    nn = 3;
    NE = zeros(1,nElem,nn);
    NE(1,:,:) = E;
    fileID = fopen('Mesh_1_26_2023.gri','w');
    fprintf(fileID,'%i %i %i\r\n',nNode,nElemTot,Dim);
    for i=1:nNode
        fprintf(fileID,'%1.4f %1.4f\r\n',N(i,:));
    end
    fprintf(fileID,'%i\r\n',nBGroup);
    for i=1:nBGroup
        fprintf(fileID,'%i %i %s\r\n',nBFace,nf(i),Title(i));
        for j=1:nf(i)-1
            fprintf(fileID,'%i ',NB(i,j,:));
        end
        fprintf(fileID,'%i\r\n',NB(i,nf(i),:));
    end
    fprintf(fileID,'%i %i %s\r\n',nElem,Order,Basis);
    for j=1:nElem
        fprintf(fileID,'%i %i %i\r\n',NE(1,j,:));
    end 
    fclose(fileID);
end

function E=triccw(E,N)
    for i=1:height(E)
        b = ccw(N(E(i,1),:),N(E(i,2),:),N(E(i,3),:));
        if b == 0
            E(i,:) = flip(E(i,:),2);
        end
    end
end
function f=ccw(a,b,c)
    f = (c(2)-a(2))*(b(1)-a(1)) > (b(2)-a(2))*(c(1)-a(1));
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
    for elem = 1:nelem,
        nv = E2N(elem,1:3);
        for edge = 1:3,
            n1 = nv(mod(edge  ,3)+1);
            n2 = nv(mod(edge+1,3)+1);
            if (H(n1,n2) == 0), % edge hit for the first time
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
