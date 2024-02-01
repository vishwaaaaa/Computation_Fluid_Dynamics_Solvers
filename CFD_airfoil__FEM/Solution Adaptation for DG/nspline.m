function N = nspline(N0,flag)
%flag = 1 for flap (left airfoil), flag = 2 for main (middle airfoil),
%flag = 3 for slat (right airfoil)
close all;
if flag == 1
    %XY = readmatrix("C:\Users\benba\Downloads\slat.txt");
    %XY = readmatrix("slat.txt");
    XY = readmatrix("flap.txt");
elseif flag == 2
    %XY = readmatrix("C:\Users\benba\Downloads\main.txt");
    XY = readmatrix("main.txt");
else
    %XY = readmatrix("C:\Users\benba\Downloads\flap.txt");
    %XY = readmatrix("flap.txt");
    XY = readmatrix("slat.txt");
end
%figure(1)
%plot(XY(:,1)',XY(:,2)','x',N0(:,1)',N0(:,2)','o')
[ppx,ppy]=spline2d(XY(:,1),XY(:,2));
syms z
for i=1:height(N0)
    [~,i1(i)] = min(sqrt((N0(i,1)-XY(:,1)).^2+(N0(i,2)-XY(:,2)).^2));
    ii = [mod(i1(i)-2,height(XY))+1 mod(i1(i),height(XY))+1];
    [~,i2(i)] = min(sqrt((N0(i,1)-XY(ii,1)).^2+(N0(i,2)-XY(ii,2)).^2));
    i2(i) = ii(i2(i));
    if i1(i) > i2(i)
        r = i1(i); i1(i) = i2(i); i2(i) = r;
    end
    ds = ppx.breaks(i2(i))-ppx.breaks(i1(i));
    ax = ppx.coefs(i1(i),:); ay = ppy.coefs(i1(i),:);
    sx = ax(1)*z^3+ax(2)*z^2+ax(3)*z+ax(4);
    sy = ay(1)*z^3+ay(2)*z^2+ay(3)*z+ay(4);
    d2 = (sx-N0(i,1))^2 + (sy-N0(i,2))^2;
    s(i) = min(abs(double(vpasolve(diff(d2)==0,z))));
    s(i) = s(i)+ppx.breaks(i1(i));
end
N(:,1) = ppval(ppx,s);
N(:,2) = ppval(ppy,s);
%figure(2)
%plot(XY(:,1)',XY(:,2)','x',N(:,1)',N(:,2)','o')     
end