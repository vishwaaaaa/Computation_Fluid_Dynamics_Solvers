function [I2E, B2E, E2N, V, In, Bn, A, Ee] = Task2_Task3_full(grifile)
% This function performs all the post-processing of the mesh. It identifies
% the interior edges and boundary edges for each node, the corresponding
% normal vectors, the area of each node, and verifies that the error on the
% element is 0
%
% INPUT : grifile str for standard .gri file
% OUTPUT: I2E = [niedge x 4] array giving (elemL, faceL, elemR, faceR)
%              information for each interior edge
%         B2E = [nbedge x 3] array giving (elem, face, bgroup)
%              information for each boundary edge
%         In = [niedge x 2] array giving (nx, ny) for each interior edge,
%              normalized
%         Bn = [nbedge x 2] array giving (nx, ny) for each boundary edge,
%              normalized
%         A = [nelem] array of the area of each triangular mesh
%         Ee = [nelem x 3] array mapping triangles to nodes

    [V, E2N, edge_info] = E2N_func(grifile);    % reading in the correct .gri file
    [Ee, norms] = verify(E2N, V);
    nelem = size(E2N,1);             % number of elements
    nnode = max(max(E2N));           % number of nodes
    H = sparse(nnode,nnode);         % Create a hash list to identify edges
    H2 = sparse(nnode,nnode);        % Create a hash list to identify faces
    H3 = sparse(nnode,nnode);
    I2E = zeros(ceil(nelem*3/2), 4); % (over) allocate interior edge array
    B2E = zeros(ceil(nelem*3/2), 3); % (over) allocate boundary edge array
    In = zeros(ceil(nelem*3/2), 2);  % (over) allocate interior edge normal array
    Bn = zeros(ceil(nelem*3/2), 2);  % (over) allocate boundary edge normal array
    A = zeros(1, nelem);             % allocate area array
    niedge = 0;                      % number of interior edges (running total)
    nbedge = 0;
    nedge = 0;
    
    % Loop over elements and identify all edges
    for elem = 1:nelem
        nv = E2N(elem,1:3);
        ledge = zeros(1,3);
        for edge = 1:3
            n1 = nv(mod(edge  ,3)+1);
            n2 = nv(mod(edge+1,3)+1);
            nedge = nedge+1;
            % calculating the length of each edge around the triangular cell
            ledge(edge) = sqrt((V(n2, 1) - V(n1, 1))^2 + (V(n2, 2) - V(n1, 2))^2);
            
            % calculating the area of a non-equilateral triangular cell
            s = (ledge(1) + ledge(2) + ledge(3))/2;                   % semi-perimeter
            A(elem) = sqrt(s*(s-ledge(1))*(s-ledge(2))*(s-ledge(3)));
            
            if (H(n1,n2) == 0) % edge hit for the first time
                % could be a boundary or interior; assume boundary
                H(n1,n2) = elem;  H(n2,n1) = elem;
                H2(n1, n2) = edge; H2(n2, n1) = edge;
                H3(n1, n2) = nedge; H3(n2, n1) = nedge;
            else % this is an interior edge, hit for the second time
                oldelem = H(n1,n2);
                oldedge = H2(n1, n2);
                if (oldelem < 0), error 'Mesh input error'; end
                niedge = niedge+1;
                elemL = min([oldelem, elem]);
                    if elemL == oldelem
                        faceL = oldedge;
                        elemR = elem;
                        faceR = edge;
                    else
                        elemL = elem; faceL = edge; elemR = elem; faceR = edge;
                    end
                I2E(niedge,:) = [elemL, faceL, elemR, faceR];
                H(n1,n2) = -1;  H(n2,n1) = -1;
                
                In(niedge,:) = norms(niedge,:);
                
            end
        end
    end
    
    I2E = I2E(1:niedge,:);  % clip I2E
    In = In(1:niedge,:);    % clip In
       
    % find boundary edges
    [I,J] = find(triu(H)>0);
    for b = 1:size(I,1) 
        nbedge = nbedge + 1;
        B2E(b, 1) = H(I(b),J(b));
        B2E(b, 2) = H2(I(b), J(b));
        face_idx = find(edge_info(:,3) == H2(I(b), J(b)));
        update = edge_info(face_idx,:);
        idx = find(update(:,4) == H(I(b),J(b)));
        B2E(b, 3) = update(idx, 5);
        
        n_edge = H3(I(b), J(b));
        Bn(nbedge, :) = norms(n_edge,:);
       
    end
    
    B2E = B2E(1:nbedge,:);  % clip B2E
    Bn = Bn(1:nbedge,:);   % clip Bn
    
end




