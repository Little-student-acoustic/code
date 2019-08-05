function U=vvin(r,theta,n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3,rp,rhop,cp)
i=sqrt(-1);
k3=omega/c3;
syms x1  x2
d2besselj=matlabFunction(diff(besselj(n,x1),x1,2));
dbesselj=matlabFunction(diff(besselj(n,x2),x2));
f1=1-(c3^(2)*rho3)/(cp^(2)*rhop);
f2=2*(rhop-rho3)/(2*rhop+rho3);
if n==0
    epsilon=1;
else
    epsilon=2;
end
u=(1/(rho3*omega))*i^(n)*epsilon*k3^(2)*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*d2besselj(k3.*r).*cos(n.*theta)*2.*((1/(rho3*omega))*i^(n)*epsilon*k3*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbesselj(k3.*r).*cos(n.*theta));
p3=i^(n)*epsilon*k3*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbesselj(k3.*r).*cos(n.*theta)*2.*(i^(n)*epsilon*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*besselj(n,k3.*r).*cos(n.*theta));
%(-n/(r*rho3*omega^(2)))*i^(n)*epsilon*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*besselj(n,k3.*r).*sin(n.*theta)+
uu=(-n/(rho3*omega))*i^(n)*epsilon*k3*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbesselj(k3.*r).*sin(n.*theta)*2.*((1/(rho3*omega))*i^(n)*epsilon*k3*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbesselj(k3.*r).*cos(n.*theta));
pp3=i^(n)*epsilon*(-n)*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*besselj(n,k3.*r).*sin(n.*theta)*2.*(i^(n)*epsilon*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*besselj(n,k3.*r).*cos(n.*theta));
% U=-2*pi*rp^(3)*rho3*(p3*(f1/3*rho3^(2)*c3^(2))-u*(f2/2))*2.*r.*sin(theta)-2*pi*rp^(3)*rho3*(pp3*(f1/3*rho3^(2)*c3^(2))-uu*(f2/2)).*(r.*cos(theta)./(r.*cos(theta).^(2)+r.*sin(theta).^(2)));
U=-2*pi*rp^(3)*rho3*(p3*(f1/(3*rho3^(2)*c3^(2)))-u*(f2/2)).*(1./sin(theta))-2*pi*rp^(3)*rho3*(pp3*(f1/(3*rho3^(2)*c3^(2)))-uu*(f2/2)).*(1./(cos(theta).*r));
end

