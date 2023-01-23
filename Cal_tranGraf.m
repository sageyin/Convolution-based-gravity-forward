function [gg]=Cal_tranGraf(x,y,z,x0,y0,z0,a,b,c,den,Style)
% Analytic formula method for calculating the gravity of a cube (subfunction)
% Editor：Xianzhe Yin 2022/9/05    
% References：Li, X., & Chouteau, M.,1998. Three-dimensional gravity modeling in all space, Surv. Geophys.,19(4), 339-368.
%% Parameters
% ===== input =====
% x y z :Coordinates of measuring network (unit:m)
% x0 y0 z0 : Center coordinates of the field source (unit:m)
% a b c : Source geometric scale (length, width and height) (unit:m)
% den : Source residual density (unit:kg/m^3)
% Style : Type of gravitational field, including the vertical component of gravity and its tensor
% ===== out =====
% gg : the vertical component of gravity or its tensor

f=6.667e-11;  % gravitational constant
xx0={x0-a/2,x0+a/2};yy0={y0-b/2,y0+b/2};zz0={z0-c/2,z0+c/2};
gg=zeros(size(x));
switch Style
    case 'gz' % the vertical component of the gravity field
        for l=1:2
            for m=1:2
                for n=1:2
                    const0=(-1)^(l)*(-1)^(m)*(-1)^(n);
                    xx=x-xx0{l};yy=y-yy0{m};zz=z-zz0{n};
                    r=sqrt(xx.^2+yy.^2+zz.^2);
                    g=-f*den*const0*(xx.*mylog(yy+r)+yy.*mylog(xx+r)+2*zz.*myatan((xx+yy+r)./zz));
                    gg=gg+g;
                end
            end
        end
    case 'gxx'
        for l=1:2
            for m=1:2
                for n=1:2
                    const0=(-1)^(l)*(-1)^(m)*(-1)^(n);
                    xx=x-xx0{l};yy=y-yy0{m};zz=z-zz0{n};
                    r=sqrt(xx.^2+yy.^2+zz.^2);
                    g=f*den*const0*(-2*myatan((yy+zz+r)./xx));
                    gg=gg+g;
                end
            end
        end
    case 'gyy'
        for l=1:2
            for m=1:2
                for n=1:2
                    const0=(-1)^(l)*(-1)^(m)*(-1)^(n);
                    xx=x-xx0{l};yy=y-yy0{m};zz=z-zz0{n};
                    r=sqrt(xx.^2+yy.^2+zz.^2);
                    g=f*den*const0*(-2*myatan((xx+zz+r)./yy));
                    gg=gg+g;
                end
            end
        end
    case 'gzz'
        for l=1:2
            for m=1:2
                for n=1:2
                    const0=(-1)^(l)*(-1)^(m)*(-1)^(n);
                    xx=x-xx0{l};yy=y-yy0{m};zz=z-zz0{n};
                    r=sqrt(xx.^2+yy.^2+zz.^2);
                    g=f*den*const0*(myatan((xx.*yy)./(zz.*r))); 
                    gg=gg+g;
                end
            end
        end
    case 'gxy'
        for l=1:2
            for m=1:2
                for n=1:2
                    const0=(-1)^(l)*(-1)^(m)*(-1)^(n);
                    xx=x-xx0{l};yy=y-yy0{m};zz=z-zz0{n};
                    r=sqrt(xx.^2+yy.^2+zz.^2);
                    g=-f*den*const0*(mylog(zz+r));
                    gg=gg+g;
                end
            end
        end
    case 'gxz'
        for l=1:2
            for m=1:2
                for n=1:2
                    const0=(-1)^(l)*(-1)^(m)*(-1)^(n);
                    xx=x-xx0{l};yy=y-yy0{m};zz=z-zz0{n};
                    r=sqrt(xx.^2+yy.^2+zz.^2);
                    g=-f*den*const0*(mylog(yy+r));
                    gg=gg+g;
                end
            end
        end
    case 'gyz'
        for l=1:2
            for m=1:2
                for n=1:2
                    const0=(-1)^(l)*(-1)^(m)*(-1)^(n);
                    xx=x-xx0{l};yy=y-yy0{m};zz=z-zz0{n};
                    r=sqrt(xx.^2+yy.^2+zz.^2);
                    g=-f*den*const0*(mylog(xx+r));
                    gg=gg+g;
                end
            end
        end
end

%% avoid singularities 
function y = mylog(x)
    y = log(x);
    y(x==0) = -10^(17);
end

function y = myatan(x)
    y = atan(x);
    y(isnan(x)) = atan2(1,0);
end





end