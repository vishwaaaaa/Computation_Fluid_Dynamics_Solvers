function [Ee, norms] = verify(E2N, V)
    
    nelem = size(E2N,1);
    norms = zeros(ceil(nelem*3/2), 2);
    Ee = zeros(nelem,2);
    niedge = 0;
    
    for elem = 1:nelem
        
        nv = E2N(elem,1:3);
        
        for edge = 1:3
            niedge = niedge + 1;
            
            n1 = nv(mod(edge  ,3)+1);
            n2 = nv(mod(edge+1,3)+1);
            
            ledge = sqrt((V(n2, 1) - V(n1, 1))^2 + (V(n2, 2) - V(n1, 2))^2);
            nx = (V(n2, 2) - V(n1, 2))/ledge;
            ny = -(V(n2, 1) - V(n1, 1))/ledge;
            
            norms(niedge, 1) = nx;
            norms(niedge, 2) = ny;
            
            Ee(elem,1) = Ee(elem,1) + nx*ledge;
            Ee(elem,2) = Ee(elem,2) + ny*ledge;
            
        end
        
    end

end