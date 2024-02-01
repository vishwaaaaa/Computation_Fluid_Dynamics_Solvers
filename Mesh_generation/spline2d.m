function [ppx,ppy]=spline2d(x,y)
    n = length(x);
    s = zeros(1,n);
    [xq,wq] = quad1();
    for i=1:n-1
        s(i+1) = s(i) + norm([x(i+1)-x(i); y(i+1)-y(i)]);
    end
    ppx = spline(s,x);
    ppy = spline(s,y);
    st = zeros(1,n);
    for i=1:n-1
        ds = s(i+1)-s(i);
        tq = 1/2*(xq+1)*ds;
        %Use with cubic
        xs = 3*ppx.coefs(i,1)*tq.^2 + 2*ppx.coefs(i,2)*tq+ ppx.coefs(i,3);
        ys = 3*ppy.coefs(i,1)*tq.^2 + 2*ppy.coefs(i,2)*tq + ppy.coefs(i,3);
        %Use with quad
        %xs = 2*ppx.coefs(i,1)*tq + ppx.coefs(i,2);
        %ys = 2*ppy.coefs(i,1)*tq + ppy.coefs(i,2);
        f = sqrt(xs.^2 + ys.^2);
        st(i+1) = st(i) + sum(f.*wq)/2*ds;
    end
    err = sum(abs(st-s))/(s(end)-s(1));
    while err > 1e-2
        s = st;
        ppx = spline(s,x);
        ppy = spline(s,y);
        st = zeros(1,n);
        for i=1:n-1
            ds = s(i+1)-s(i);
            tq = 1/2*(xq+1)*ds;
            %Use with cubic
            xs = 3*ppx.coefs(i,1)*tq.^2 + 2*ppx.coefs(i,2)*tq+ ppx.coefs(i,3);
            ys = 3*ppy.coefs(i,1)*tq.^2 + 2*ppy.coefs(i,2)*tq + ppy.coefs(i,3);
            %Use with quad
            %xs = 2*ppx.coefs(i,1)*tq + ppx.coefs(i,2);
            %ys = 2*ppy.coefs(i,1)*tq + ppy.coefs(i,2);
            f = sqrt(xs.^2 + ys.^2);
            st(i+1) = st(i) + sum(f.*wq)/2*ds;
        end
        err = sum(abs(st-s))/(s(end)-s(1));
    end
end
%Bonus function, does quadratic splines instead of cubic splines
function pp=qspline(s,x)
    n = length(s);
    m = zeros(1,n);
    a2(1) = 0; %zero second derivative at first point
    a0(1) = x(1);
    a1(1) = (x(2)-x(1))/(s(2)-s(1));
    m(1) = a1(1);
    m(2) = m(1);
    for i=2:n-1
        ds(i) = s(i+1)-s(i);
        a0(i) = x(i);
        a1(i) = m(i);
        a2(i) = x(i+1) - a0(i) - a1(i)*ds(i);
        m(i+1) = (a1(i) + 2*a2(i)*ds(i)); %march along arc by matching slopes
    end
    pp.form = 'pp';
    pp.breaks = s;
    pp.coefs = [a2' a1' a0'];
    pp.pieces = length(s)-1;
    pp.order = 3;
    pp.dim = 1;
end
function [x, w] = quad1()
    x = zeros(10,1);
    w = zeros(10,1);
    x(0+1) = -0.9739065285171717200779640;
    x(1+1) = -0.8650633666889845107320967;
    x(2+1) = -0.6794095682990244062343274;
    x(3+1) = -0.4333953941292471907992659;
    x(4+1) = -0.1488743389816312108848260;
    x(5+1) =  0.1488743389816312108848260;
    x(6+1) =  0.4333953941292471907992659;
    x(7+1) =  0.6794095682990244062343274;
    x(8+1) =  0.8650633666889845107320967;
    x(9+1) =  0.9739065285171717200779640;
    w(0+1) =  0.0666713443086881375935688;
    w(1+1) =  0.1494513491505805931457763;
    w(2+1) =  0.2190863625159820439955349;
    w(3+1) =  0.2692667193099963550912269;
    w(4+1) =  0.2955242247147528701738930;
    w(5+1) =  0.2955242247147528701738930;
    w(6+1) =  0.2692667193099963550912269;
    w(7+1) =  0.2190863625159820439955349;
    w(8+1) =  0.1494513491505805931457763;
    w(9+1) =  0.0666713443086881375935688;
end