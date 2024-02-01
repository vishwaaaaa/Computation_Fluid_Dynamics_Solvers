function phi = Basis(xy, p) 
x = xy(1);
y = xy(2);
rank = (p+1)*(p+2)/2;
phi = zeros(rank,1);
switch p
    case 0
        phi(1) = 1;
        
    case 1
        phi(1) = 1-x-y;
        phi(2) = x;
        phi(3) = y;
        
    case 2
        phi(1) = 1.0-3.0*x-3.0*y+2.0*x*x+4.0*x*y+2.0*y*y;
        phi(3) = -x+2.0*x*x;
        phi(6) = -y+2.0*y*y;
        phi(5) = 4.0*x*y;
        phi(4) = 4.0*y-4.0*x*y-4.0*y*y;
        phi(2) = 4.0*x-4.0*x*x-4.0*x*y;
        
    case 3
        phi(1) = 1.0-11.0/2.0*x-11.0/2.0*y+9.0*x*x+18.0*x*y+9.0*y*y-9.0/2.0*x*x*x-27.0/2.0*x*x*y-27.0/2.0*x*y*y-9.0/2.0*y*y*y;
        phi(4)= x-9.0/2.0*x*x+9.0/2.0*x*x*x;
        phi(10) = y-9.0/2.0*y*y+9.0/2.0*y*y*y;
        phi(7) = -9.0/2.0*x*y+27.0/2.0*x*x*y;
        phi(9) = -9.0/2.0*x*y+27.0/2.0*x*y*y;
        phi(8) = -9.0/2.0*y+9.0/2.0*x*y+18.0*y*y-27.0/2.0*x*y*y-27.0/2.0*y*y*y;
        phi(5) = 9.0*y-45.0/2.0*x*y-45.0/2.0*y*y+27.0/2.0*x*x*y+27.0*x*y*y+27.0/2.0*y*y*y;
        phi(2) = 9.0*x-45.0/2.0*x*x-45.0/2.0*x*y+27.0/2.0*x*x*x+27.0*x*x*y+27.0/2.0*x*y*y;
        phi(3) = -9.0/2.0*x+18.0*x*x+9.0/2.0*x*y-27.0/2.0*x*x*x-27.0/2.0*x*x*y;
        phi(6) = 27.0*x*y-27.0*x*x*y-27.0*x*y*y;
        
    otherwise
        cout << "Invalid Order \n";
        phi(0,0) = 0;
        
end
end
