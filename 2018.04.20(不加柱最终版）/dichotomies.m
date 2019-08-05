%¶þ·Ö·¨
function omega=dichotomies(omega1,omega2,n,rho1,rho2,rho3,a,b,c1,cl2,ct2,c3,precision)
while omega2-omega1>precision
    y1=coefficient(n,rho1,rho2,rho3,a,b,omega1,c1,cl2,ct2,c3);
    y2=coefficient(n,rho1,rho2,rho3,a,b,omega2,c1,cl2,ct2,c3);
    omega0=(omega1+omega2)/2;
    y0=coefficient(n,rho1,rho2,rho3,a,b,omega0,c1,cl2,ct2,c3);
    if y1*y0<0
        omega2=omega0;
    elseif y2*y0<0
        omega1=omega0;
    elseif y0==0
        break
    elseif y1==0
        omega2=omega1;
        break
    elseif y2==0
        omega1=omega2;
        break
    end
end
omega=(omega1+omega2)/2;
end

