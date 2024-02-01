%%%%% AE623 Plotting Code %%%%%

close all
clear
clc

% nodes = load('UnitCurveN2.txt');
% elem = load('UnitCurveE2.txt');

% nodes = load('NewUnitMeshN.txt');
% elem = load('NewUnitMeshE.txt');
% B2E = load('NewUnitMeshB2E.txt');
% I2E = load('NewUnitMeshI2E.txt');

% nodes = load('c0_q1_nodes.txt');
% elem = load('c0_q1_E2N.txt');
% B2E = load('c0_q1_B2E.txt');
% I2E = load("c0_q1_I2E.txt");
% 
% states = load('c0_p1_q1_U.txt');


nodes = load('c0_q2_nodes.txt');
elem = load('c0_q22_E2N.txt');
B2E = load('c0_q1_B2E.txt');
I2E = load("c0_q1_I2E.txt");

states = load('c0_p2_q2_U.txt');



%%%%%% Manually Switch for Solution Input %%%%
p = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%seperate out cells with boundaries, and which edge is the boundary edge
wallIndex = 1;
flapIndex = 1;
mainIndex = 1;
slatIndex = 1;
for i = 1:length(B2E)
    if B2E(i,3) == -1
        wallCells(wallIndex,1) = B2E(i,1);
        wallCells(wallIndex,2) = B2E(i,2);
        wallIndex = wallIndex+1;
    elseif B2E(i,3) == -2
        flapCells(flapIndex,1) = B2E(i,1);
        flapCells(flapIndex,2) = B2E(i,2);
        flapIndex = flapIndex+1;
    elseif B2E(i,3) == -3
        mainCells(mainIndex,1) = B2E(i,1);
        mainCells(mainIndex,2) = B2E(i,2);
        mainIndex = mainIndex+1;
    elseif B2E(i,3) == -4
        slatCells(slatIndex,1) = B2E(i,1);
        slatCells(slatIndex,2) = B2E(i,2);
        slatIndex = slatIndex+1;
    end
end

%combine and sort by index low to high, that way we dont need to do
%complicated "find" functions when looping through the mesh to find the
%boundary cells
boundaryEdges = [flapCells; slatCells; mainCells];
boundaryEdges = sortrows(boundaryEdges);
boundaryFlag = 1;


%generate free stream flow values
Minf = 0.25;
alpha = 8;
gamma = 1.4;

rho0 = 1;
c0 = 1;
p0 = c0^2*rho0/gamma;


rhoInf = rho0 * (1+(gamma-1)/2 * Minf^2)^(-1/(gamma-1));
pInf = p0 * (1+(gamma-1)/2 * Minf^2)^(-gamma/(gamma-1));
aInf = sqrt(gamma * pInf/rhoInf);

uInf = Minf*aInf*cosd(alpha);
vInf = Minf*aInf*sind(alpha);
velInf = [uInf vInf];

rhoEInf = pInf/(gamma-1) + 1/2*rhoInf*(uInf^2+vInf^2);

chord = 1;

qinf = 1/2*rhoInf*norm(velInf)^2*chord;




%read in solution values for varying p/q values
rank = (p+1)*(p+2)/2;
edges = length(elem);
Mach = zeros(1,size(states,1));

%read in states
rho_unrolled = states(:,1);
rhou_unrolled = states(:,2);
rhov_unrolled = states(:,3);
rhoE_unrolled = states(:,4);


%want to form a list of states per each element @ lagrange points
rho = zeros(length(rho_unrolled)/rank,rank);
rhou = zeros(length(rhou_unrolled)/rank,rank);
rhov = zeros(length(rhov_unrolled)/rank,rank);
rhoE = zeros(length(rhoE_unrolled)/rank,rank);

%step through unrolled state and assign all points per element to one row
for i = 1:size(rho,1)
    for j = 1:rank
        iloc = (i-1)*rank;
        rho(i,j) = rho_unrolled(iloc+j);
        rhou(i,j) = rhou_unrolled(iloc+j);
        rhov(i,j) = rhov_unrolled(iloc+j);
        rhoE(i,j) = rhoE_unrolled(iloc+j);
    end
end
% 
% hold on
% for index = 1:size(elem,1)
%     line([nodes(elem(index,1),1) nodes(elem(index,2),1)] , [nodes(elem(index,1),2) nodes(elem(index,2),2)])
%     line([nodes(elem(index,1),1) nodes(elem(index,3),1)] , [nodes(elem(index,1),2) nodes(elem(index,3),2)])
%     line([nodes(elem(index,2),1) nodes(elem(index,3),1)] , [nodes(elem(index,2),2) nodes(elem(index,3),2)])
% end




Drag = 0;
Lift = 0;

 %%%%%%%%%%%
    N = 5; % number of intervals for plotting
 %%%%%%%%%%%

%pre allocate a 3d matrix of solutions @ plotting points in ref space 
% so we can plot after the main loop, speeds up code dramatically
 V = zeros((N+1)*(N+2)/2, 2, size(elem,1));
    
%step through every cell in the mesh
for elementIndex = 1:size(elem,1)
    elementIndex
    element = elem(elementIndex,:);

    xi = linspace(0,1,N+1);
    eta = linspace(0,1,N+1);
   
    M = zeros(N+1,N+1); % mapping matrix
    k = 0;

    for j = 0:N
        for i = 0:(N-j)
            k = k + 1;
            V(k,1,elementIndex) = xi(i+1);
            V(k,2,elementIndex) = eta(j+1);
            M(i+1,j+1) = k; % node number in spot i,j
        end
    end

    %compute basis function solutions @ plotting points in ref space
    for i = 1:size(V,1)
        phiMat(i,:) = TriLagrange(V(i,1,elementIndex),V(i,2,elementIndex),p);
    end

    % evaluate the plotting point state interpolation
    rhoPlot(elementIndex,:) = phiMat*rho(elementIndex,:)';
    rhouPlot(elementIndex,:) = phiMat*rhou(elementIndex,:)';
    rhovPlot(elementIndex,:) = phiMat*rhov(elementIndex,:)';
    rhoEPlot(elementIndex,:) = phiMat*rhoE(elementIndex,:)';
    
    %solve for state variables at each plotting point

    for i = 1:size(rhoPlot(elementIndex,:),2)
        uplot(elementIndex,i) = rhouPlot(elementIndex,i)./rhoPlot(elementIndex,i);
        vplot(elementIndex,i) = rhovPlot(elementIndex,i)./rhoPlot(elementIndex,i);

        pressure(elementIndex,i) = (gamma-1)*(rhoEPlot(elementIndex,i)...
            - 1/2 * rhoPlot(elementIndex,i)*...
            (uplot(elementIndex,i)^2 + vplot(elementIndex,i)^2));

        c(elementIndex,i) = sqrt(gamma * pressure(elementIndex,i)/rhoPlot(elementIndex,i));
        mach(elementIndex,i) = norm([uplot(elementIndex,i) vplot(elementIndex,i)])/c(elementIndex,i);
    end
      
      
    
    

    % E = subelement matrix
    E(:,:,elementIndex) = zeros(N*N,3);
    k = 0;
    for j = 1:N
        for i = 1:(N+1-j)
            k = k+1;
            E(k,:,elementIndex) = [M(i,j), M(i+1,j), M(i,j+1)];
            if (i < N+1-j)
                k = k+1;
                E(k,:,elementIndex) = [M(i+1,j), M(i+1,j+1), M(i,j+1)];
            end
        end
    end

    % add a V mapping for curved elements

    %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = q; % geometry order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%


    for j=1:length(element)
    xyQ(j,:) = nodes(element(j),:);
    end
    

    V(:,:,elementIndex) = Ref2Glob(V(:,:,elementIndex), Q, xyQ);
   
for edge = 1:3
    if edge == 1
        xi_sigma = -1; eta_sigma = 1;
    elseif edge == 2
        xi_sigma = 0; eta_sigma = -1;
    elseif edge == 3
        xi_sigma = 1; eta_sigma = 0;
    end
    [xq1, wq1] = quad1d(p,q);
    xref = RefEdge2RefElem(edge, xq1);
    lengths = 0;
    for i=1:size(xref,1)
       J = ElemJacobian(xref(i,:), Q, xyQ);
       stan = J(:,1)*xi_sigma + J(:,2)*eta_sigma; % tangent * ds/dsigma
       nvec = [stan(2), -stan(1)];
       lengths = lengths+norm(nvec);
    end
lengthOut(elementIndex,edge) = lengths/i;
end



    % for our presorted edge flags we found in the file read in portion
    if elementIndex == boundaryEdges(boundaryFlag,1)
        edge = boundaryEdges(boundaryFlag,2);
            if edge == 1
                xi_sigma = -1; eta_sigma = 1;
            elseif edge == 2
                xi_sigma = 0; eta_sigma = -1;
            elseif edge == 3
                xi_sigma = 1; eta_sigma = 0;
            end
            [xq1, wq1] = quad1d(p,q);
            xref = RefEdge2RefElem(edge, xq1);
            
               %compute basis function solutions @ edge quad points in ref space
                    for i = 1:size(xref,1)
                        phiMatEdge(i,:) = TriLagrange(xref(i,1),xref(i,2),p);
                    end
                
                    % evaluate the edge quad point state interpolation
                    rhoEdge = phiMatEdge*rho(elementIndex,:)';
                    rhouEdge = phiMatEdge*rhou(elementIndex,:)';
                    rhovEdge = phiMatEdge*rhov(elementIndex,:)';
                    rhoEEdge = phiMatEdge*rhoE(elementIndex,:)';

                    for i = 1:length(rhoEdge)
                        uEdge(i) = rhouEdge(i)/rhoEdge(i);
                        vEdge(i) = rhovEdge(i)/rhoEdge(i);
                        vel(i,:) = [uEdge(i) vEdge(i)];
                       
                        pEdge(i) = (gamma-1)*(rhoEEdge(i) - 1/2 * rhoEdge(i)*(uEdge(i)^2 + vEdge(i)^2));
                        cEdge(i) = sqrt(gamma * pEdge(i)/rhoEdge(i));
                    end
              
            lengths = 0;
            totX = 0;
            totY = 0;
            for i=1:size(xref,1)
               J = ElemJacobian(xref(i,:), Q, xyQ);
               stan = J(:,1)*xi_sigma + J(:,2)*eta_sigma; % tangent * ds/dsigma
               nvec = [stan(2), -stan(1)];
               nhat = nvec./norm(nvec);

               normalVelx = vel(i,1) - dot(nvec,vel(i,:))*nvec(1);
               normalVely = vel(i,2) - dot(nvec,vel(i,:))*nvec(2);

               pressX = nhat(1)*pEdge(i)*wq1(i);
               pressY = nhat(2)*pEdge(i)*wq1(i);
                
               totX = totX + pressX;
               totY = totY + pressY;

               normalMach(boundaryFlag,i) = norm([normalVelx normalVely])/cEdge(i);
               wallPressure(i) = pEdge(i)*wq1(i);
               
    
               lengths = lengths+norm(nvec);
            end
            
               Drag = Drag + totX*cosd(alpha)*lengths;
               Lift = Lift + totY*sind(alpha)*lengths;
             
              coeffPress(boundaryFlag) = (sum(wallPressure) - pInf)/qinf;

              boundaryFlag = boundaryFlag+1;
                  if boundaryFlag > length(boundaryEdges)
                   boundaryFlag = 1;    
                    %once we reach max boundaries, reset to one to just ignore and not
                    %array out of bounds upon next loop
                    end
    end
end


%%

% figure(1)
% hold on
% axis([-0.3 1.4 -0.6 0.6]);
% for el = 1:size(elem,1)
%     el
%     for k = 1:size(E,1)
%         X = V(E(k,:,el), :, el);
%         I = [1,2,3];
%         patch(X(I,1), X(I,2), [(mach(el,E(k,:,el)))]','LineStyle','none'); hold on; 
%     end
% end
% 
% %axis equal
% colormap('jet')
% colorbar
% xlabel('x-Coordinate')
% ylabel('y-Coordinate')
% title("Mach Contours")
% subtitle("Coarse Mesh, p = " + p + ", q = " + q)

%%

for i = 1:length(coeffPress)
xvector(i) =  nodes(elem(boundaryEdges(i,1),boundaryEdges(i,2)),1);
end
figure(2)
hold on
cpPlotting = [xvector; (coeffPress)];
cpPlotFinal = sortrows(cpPlotting');

posIndex = 1;
negIndex = 1;
for i = 1:length(cpPlotFinal)
    if cpPlotFinal(i,2) > 0
        cpPos(posIndex,:) = cpPlotFinal(i,:);
        posIndex = posIndex+1;
    else
        cpNeg(negIndex,:) = cpPlotFinal(i,:);
        negIndex = negIndex+1;
    end
end

for i = 2:size(cpPos,1)
    line([cpPos(i-1,1) cpPos(i,1)],[cpPos(i,2) cpPos(i,2)],'LineWidth',3)
end
for i = 2:size(cpNeg,1)
    line([cpNeg(i-1,1) cpNeg(i,1)],[cpNeg(i,2) cpNeg(i,2)],'LineWidth',3)
end
set(gca, 'YDir','reverse')
xlabel('x-Coordinate')
ylabel('Coefficient of Pressure')
title('Coefficient of Pressure vs x-axis')
subtitle("Coarse Mesh, p = " + p + ", q = " + q)
grid on

Lift
abs(Drag) 

CL = Lift/(qinf)
CD = abs(Drag)/(qinf)

%%

for i = 1:size(mach,1)
    avgMach(i) = mean(mach(i,:));
end
for i = 1:size(normalMach,1)
    mapToElementBoundaryMach = boundaryEdges(i,1);
    avgBoundaryMach = mean(normalMach(i,:));
    BMout(i,:) = [mapToElementBoundaryMach,avgBoundaryMach];
end


r = -0.5;
for i=1:size(I2E,1)
    cellL = I2E(i,1);
    cellR = I2E(i,3);
    edgeL = I2E(i,2);
    edgeR = I2E(i,4);

    lengthValue = lengthOut(cellL,edgeL);
    machDiff = abs(avgMach(cellL) - avgMach(cellR));
    interiorE = lengthValue^r * machDiff;

    interiorEout(i,:) = [cellL cellR interiorE];
end

for i=1:size(B2E,1)
     boundaryCellL = B2E(i,1);
     boundaryR = B2E(i,3);
     boundaryEdgeL = B2E(i,2);
     lengthValue = lengthOut(boundaryCellL,boundaryEdgeL);

     if boundaryR == -1
         boundaryE = 0;
         boundaryEout(i,:) = [boundaryCellL boundaryR boundaryE];
     else
         BMindex = find(BMout(:,1) == boundaryCellL);
         boundaryE = lengthValue^r * (BMout(BMindex,2));
         boundaryEout(i,:) = [boundaryCellL boundaryR boundaryE];
     end

     
end

adaptiveOutput = [interiorEout; boundaryEout];
%save('adaptMetric_q1_p0.txt', 'adaptiveOutput' ,'-ascii')





%--------------------------------
function xglob = Ref2Glob(xref, Q, xyQ)
xglob = xref;
for n = 1:size(xref,1)
phi = TriLagrange(xref(n,1), xref(n,2), Q);
xy = [0,0];
for i = 1:size(xyQ,1)
xy = xy + phi(i)*xyQ(i,:);
end
xglob(n,:) = xy;
end
end
%--------------------------------
function xref = RefEdge2RefElem(edge, xedge)
sigma = zeros(length(xedge), 1);
sigma(:,1) = xedge; Z = zeros(size(sigma));
if (edge == 1)
xref = [1-sigma, sigma];
elseif (edge == 2)
xref = [Z, 1-sigma];
elseif (edge == 3)
xref = [sigma, Z];
else
error('edge out of bounds');
end
end
%--------------------------------
function J = ElemJacobian(xref, Q, xyQ)
[phix, phiy] = GradTriLagrange(xref(1),xref(2),Q);
J = zeros(2,2);
for i = 1:size(xyQ,1)
J = J + [xyQ(i,:)'*phix(i), xyQ(i,:)'*phiy(i)];
end
end
%--------------------------------
function phi = TriLagrange(x,y,p) 
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
        phi(1) = 0;
        
end
end
%-------------------------------
function [gphix, gphiy] = GradTriLagrange(x,y,p)
switch p
    case 0
        gphix =  0.0;
        gphiy =  0.0; 

    case 1
        gphix(1) =  -1.0;
        gphix(2) =  1.0;
        gphix(3) =  0.0;
        gphiy(1) =  -1.0;
        gphiy(2) =  0.0;
        gphiy(3) =  1.0;

    case 2
        gphix(1) =  -3.0+4.0*x+4.0*y;
        gphix(3) =  -1.0+4.0*x;
        gphix(6) =  0.0;
        gphix(5) =  4.0*y;
        gphix(4) =  -4.0*y;
        gphix(2) =  4.0-8.0*x-4.0*y;
        gphiy(1) =  -3.0+4.0*x+4.0*y;
        gphiy(3) =  0.0;
        gphiy(6) =  -1.0+4.0*y;
        gphiy(5) =  4.0*x;
        gphiy(4) =  4.0-4.0*x-8.0*y;
        gphiy(2) =  -4.0*x;

    case 3
        gphix(1) =  -11.0/2.0+18.0*x+18.0*y-27.0/2.0*x*x-27.0*x*y-27.0/2.0*y*y;
        gphix(4) =  1.0-9.0*x+27.0/2.0*x*x;
        gphix(10) =  0.0;
        gphix(7) =  -9.0/2.0*y+27.0*x*y;
        gphix(9) =  -9.0/2.0*y+27.0/2.0*y*y;
        gphix(8) =  9.0/2.0*y-27.0/2.0*y*y;
        gphix(5) =  -45.0/2.0*y+27.0*x*y+27.0*y*y;
        gphix(2) =  9.0-45.0*x-45.0/2.0*y+81.0/2.0*x*x+54.0*x*y+27.0/2.0*y*y;
        gphix(3) =  -9.0/2.0+36.0*x+9.0/2.0*y-81.0/2.0*x*x-27.0*x*y;
        gphix(6) =  27.0*y-54.0*x*y-27.0*y*y;
        
        gphiy(1) =  -11.0/2.0+18.0*x+18.0*y-27.0/2.0*x*x-27.0*x*y-27.0/2.0*y*y;
        gphiy(4) =  0.0;
        gphiy(10) =  1.0-9.0*y+27.0/2.0*y*y;
        gphiy(7) =  -9.0/2.0*x+27.0/2.0*x*x;
        gphiy(9) =  -9.0/2.0*x+27.0*x*y;
        gphiy(8) =  -9.0/2.0+9.0/2.0*x+36.0*y-27.0*x*y-81.0/2.0*y*y;
        gphiy(5) =  9.0-45.0/2.0*x-45.0*y+27.0/2.0*x*x+54.0*x*y+81.0/2.0*y*y;
        gphiy(2) =  -45.0/2.0*x+27.0*x*x+27.0*x*y;
        gphiy(3) =  9.0/2.0*x-27.0/2.0*x*x;
        gphiy(6) =  27.0*x-27.0*x*x-54.0*x*y;
     
    otherwise
        cout << "Invalid Order \n";
        gphix(1) = 0;
        gphiy(1) = 0;

end
end
%-------------------------------
function [xp, wp] = quad1d(p, q)
    order = 2*p+1+q-1;
     if (mod(order,2)) == 0
         order = order+1;
     end

    switch order
    case 1
        x = [0.5];
        w = [1];
    case 3
        x = [0.211324865405187, 0.788675134594813];
        w = [0.5, 0.5];
    case 5
        x = [0.112701665379258, 0.500000000000000, 0.887298334620742];
        w = [0.277777777777778, 0.444444444444444, 0.277777777777778];
    case 7
        x = [0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026];
        w = [0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727];
    case 9
        x = [0.046910077030668, 0.230765344947158, 0.500000000000000, 0.769234655052841, 0.953089922969332];
        w = [0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095];
    case 11
        x = [0.033765242898424, 0.169395306766868, 0.380690406958402, 0.619309593041598, 0.830604693233132, 0.966234757101576];
        w = [0.085662246189585, 0.180380786524069, 0.233956967286345, 0.233956967286345, 0.180380786524069, 0.085662246189585];
    case 13
        x = [0.025446043828621, 0.129234407200303, 0.297077424311301, 0.500000000000000, 0.702922575688699, 0.870765592799697, 0.974553956171379];
        w = [0.064742483084435, 0.139852695744638, 0.190915025252560, 0.208979591836735, 0.190915025252560, 0.139852695744638, 0.064742483084435];
    case 15
        x = [0.019855071751232, 0.101666761293187, 0.237233795041836, 0.408282678752175, 0.591717321247825, 0.762766204958164, 0.898333238706813, 0.980144928248768];
        w = [0.050614268145188, 0.111190517226687, 0.156853322938944, 0.181341891689181, 0.181341891689181, 0.156853322938944, 0.111190517226687, 0.050614268145188];
    case 17
        x = [0.015919880246187, 0.081984446336682, 0.193314283649705, 0.337873288298096, 0.500000000000000, 0.662126711701905, 0.806685716350295, 0.918015553663318, 0.984080119753813];
        w = [0.040637194180787, 0.090324080347429, 0.130305348201468, 0.156173538520001, 0.165119677500630, 0.156173538520001, 0.130305348201468, 0.090324080347429, 0.040637194180787];
    case 19
        x = [0.013046735741414, 0.067468316655508, 0.160295215850488, 0.283302302935376, 0.425562830509184, 0.574437169490816, 0.716697697064624, 0.839704784149512, 0.932531683344492, 0.986953264258586];
        w = [0.033335672154344, 0.074725674575290, 0.109543181257991, 0.134633359654998, 0.147762112357376, 0.147762112357376, 0.134633359654998, 0.109543181257991, 0.074725674575290, 0.033335672154344];
    case 21
        x = [0.010885670926972, 0.056468700115952, 0.134923997212975, 0.240451935396594, 0.365228422023828, 0.500000000000000, 0.634771577976172, 0.759548064603406, 0.865076002787025, 0.943531299884048, 0.989114329073028];
        w = [0.027834283558087, 0.062790184732452, 0.093145105463867, 0.116596882295995, 0.131402272255123, 0.136462543388950, 0.131402272255123, 0.116596882295995, 0.093145105463867, 0.062790184732452, 0.027834283558087];
    case 23
        x = [0.009219682876640, 0.047941371814763, 0.115048662902848, 0.206341022856691, 0.316084250500910, 0.437383295744266, 0.562616704255734, 0.683915749499090, 0.793658977143309, 0.884951337097152, 0.952058628185237, 0.990780317123360];
        w = [0.023587668193256, 0.053469662997659, 0.080039164271673, 0.101583713361533, 0.116746268269177, 0.124573522906701, 0.124573522906701, 0.116746268269177, 0.101583713361533, 0.080039164271673, 0.053469662997659, 0.023587668193256];
    case 25
        x = [0.007908472640706, 0.041200800388511, 0.099210954633345, 0.178825330279830, 0.275753624481777, 0.384770842022433, 0.500000000000000, 0.615229157977567, 0.724246375518223, 0.821174669720170, 0.900789045366655, 0.958799199611489, 0.992091527359294];
        w = [0.020242002382658, 0.046060749918864, 0.069436755109894, 0.089072990380973, 0.103908023768444, 0.113141590131449, 0.116275776615437, 0.113141590131449, 0.103908023768444, 0.089072990380973, 0.069436755109894, 0.046060749918864, 0.020242002382658];
    case 27
        x = [0.006858095651594, 0.035782558168213, 0.086399342465117, 0.156353547594157, 0.242375681820923, 0.340443815536055, 0.445972525646328, 0.554027474353672, 0.659556184463945, 0.757624318179077, 0.843646452405843, 0.913600657534883, 0.964217441831787, 0.993141904348406];
        w = [0.017559730165876, 0.040079043579880, 0.060759285343952, 0.078601583579097, 0.092769198738969, 0.102599231860648, 0.107631926731579, 0.107631926731579, 0.102599231860648, 0.092769198738969, 0.078601583579097, 0.060759285343952, 0.040079043579880, 0.017559730165876];
    case 29
        x = [0.006003740989757, 0.031363303799647, 0.075896708294786, 0.137791134319915, 0.214513913695731, 0.302924326461218, 0.399402953001283, 0.500000000000000, 0.600597046998717, 0.697075673538782, 0.785486086304269, 0.862208865680085, 0.924103291705214, 0.968636696200353, 0.993996259010243];
        w = [0.015376620998059, 0.035183023744054, 0.053579610233586, 0.069785338963077, 0.083134602908497, 0.093080500007781, 0.099215742663556, 0.101289120962781, 0.099215742663556, 0.093080500007781, 0.083134602908497, 0.069785338963077, 0.053579610233586, 0.035183023744054, 0.015376620998059];
    case 31
        x = [0.005299532504175, 0.027712488463384, 0.067184398806084, 0.122297795822498, 0.191061877798678, 0.270991611171386, 0.359198224610371, 0.452493745081181, 0.547506254918819, 0.640801775389629, 0.729008388828614, 0.808938122201322, 0.877702204177502, 0.932815601193916, 0.972287511536616, 0.994700467495825];
        w = [0.013576229705877, 0.031126761969324, 0.047579255841246, 0.062314485627767, 0.074797994408288, 0.084578259697501, 0.091301707522462, 0.094725305227534, 0.094725305227534, 0.091301707522462, 0.084578259697501, 0.074797994408288, 0.062314485627767, 0.047579255841246, 0.031126761969324, 0.013576229705877];
    case 33
        x = [0.004712262342791, 0.024662239115616, 0.059880423136507, 0.109242998051599, 0.171164420391655, 0.243654731456762, 0.324384118273062, 0.410757909252076, 0.500000000000000, 0.589242090747924, 0.675615881726938, 0.756345268543238, 0.828835579608345, 0.890757001948401, 0.940119576863493, 0.975337760884384, 0.995287737657209];
        w = [0.012074151434274, 0.027729764686994, 0.042518074158590, 0.055941923596702, 0.067568184234263, 0.077022880538405, 0.084002051078225, 0.088281352683496, 0.089723235178103, 0.088281352683496, 0.084002051078225, 0.077022880538405, 0.067568184234263, 0.055941923596702, 0.042518074158590, 0.027729764686994, 0.012074151434274];
    case 35
        x = [0.004217415789535, 0.022088025214301, 0.053698766751222, 0.098147520513738, 0.154156478469823, 0.220114584463026, 0.294124419268579, 0.374056887154247, 0.457612493479132, 0.542387506520868, 0.625943112845753, 0.705875580731421, 0.779885415536974, 0.845843521530177, 0.901852479486262, 0.946301233248778, 0.977911974785699, 0.995782584210466];
        w = [0.010808006763242, 0.024857274447485, 0.038212865127445, 0.050471022053144, 0.061277603355739, 0.070321457335325, 0.077342337563133, 0.082138241872916, 0.084571191481572, 0.084571191481572, 0.082138241872916, 0.077342337563133, 0.070321457335325, 0.061277603355739, 0.050471022053144, 0.038212865127445, 0.024857274447485, 0.010808006763242];
    case 37 
        x = [0.003796578078208, 0.019895923932585, 0.048422048192591, 0.088642671731429, 0.139516911332385, 0.199727347669160, 0.267714629312020, 0.341717950018185, 0.419820677179887, 0.500000000000000, 0.580179322820113, 0.658282049981815, 0.732285370687980, 0.800272652330841, 0.860483088667615, 0.911357328268571, 0.951577951807409, 0.980104076067415, 0.996203421921792];
        w = [0.009730894114863, 0.022407113382850, 0.034522271368821, 0.045745010811225, 0.055783322773667, 0.064376981269668, 0.071303351086803, 0.076383021032930, 0.079484421696977, 0.080527224924392, 0.079484421696977, 0.076383021032930, 0.071303351086803, 0.064376981269668, 0.055783322773667, 0.045745010811225, 0.034522271368821, 0.022407113382850, 0.009730894114863];
    case 39 
        x = [0.003435700407453, 0.018014036361043, 0.043882785874337, 0.080441514088891, 0.126834046769925, 0.181973159636742, 0.244566499024586, 0.313146955642290, 0.386107074429177, 0.461736739433251, 0.538263260566749, 0.613892925570823, 0.686853044357710, 0.755433500975414, 0.818026840363258, 0.873165953230075, 0.919558485911109, 0.956117214125663, 0.981985963638957, 0.996564299592547];
        w = [0.008807003569576, 0.020300714900193, 0.031336024167055, 0.041638370788352, 0.050965059908620, 0.059097265980759, 0.065844319224588, 0.071048054659191, 0.074586493236302, 0.076376693565363, 0.076376693565363, 0.074586493236302, 0.071048054659191, 0.065844319224588, 0.059097265980759, 0.050965059908620, 0.041638370788352, 0.031336024167055, 0.020300714900193, 0.008807003569576]; 
    end
    xp = x;
    wp = w;
end       
        
        
 