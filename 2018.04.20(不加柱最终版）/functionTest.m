%Ñ°¸ù
function [omega0,num]=functionTest(omega_min,omega_max,rho1,rho2,rho3,a,b,c1,cl2,ct2,c3,n,mode,scale)
num=0;
delta=2*pi*(1.5*10^-4)*c1/scale;
omega=omega_min:delta:omega_max;
omega1=omega(1:end-1);
omega2=omega(2:end);
omega0=zeros(mode,2);
y1=coefficient(n,rho1,rho2,rho3,a,b,omega1,c1,cl2,ct2,c3);
y2=coefficient(n,rho1,rho2,rho3,a,b,omega2,c1,cl2,ct2,c3);
y=y1.*y2;
for i=1:length(y)
    if y(i)<0 && num<mode
        num=num+1;
        omega0(num,1)=omega1(i);
        omega0(num,2)=omega2(i);
        if num==mode
            break
        end
    end
end
end

