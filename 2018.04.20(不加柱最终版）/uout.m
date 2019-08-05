function U=uout(r,theta,n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3,rp,rhop,cp)
i=sqrt(-1);
k1=omega/c1;
syms x1 x2 x3 x4
dbesselj=matlabFunction(diff(besselj(n,x1),x1));
dbesselh=matlabFunction(diff(besselj(n,x3)+i*bessely(n,x3),x3));
d2besselj=matlabFunction(diff(besselj(n,x2),x2,2));
d2besselh=matlabFunction(diff(besselj(n,x4)+i*bessely(n,x4),x4,2));
f1=1-(c3^(2)*rho3)/(cp^(2)*rhop);
f2=2*(rhop-rho1)/(2*rhop+rho1);
if n==0
    epsilon=1;
else
    epsilon=2;
end
u=(1/(rho1*omega^2))*i^(n)*epsilon*k1^(2)*cos(n.*theta).*(d2besselj(k1.*r)+(Bn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)...
    /Dn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3))*d2besselh(k1.*r))*2.*((1/(rho1*omega^2))*i^(n)*epsilon*k1*cos(n.*theta).*(dbesselj(k1.*r)+(Bn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)...
    /Dn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3))*dbesselh(k1.*r)));
p1=i^(n)*epsilon*k1*(dbesselj(k1.*r)+(Bn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)...
    /Dn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3))*dbesselh(k1.*r)).*cos(n.*theta)*2.*(i^(n)*epsilon*(besselj(n,k1.*r)+(Bn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)...
    /Dn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3))*besselh(n,1,k1.*r)).*cos(n.*theta));
%(-n/(r*rho1*omega^(2)))*i^(n)*epsilon*(besselj(n,k1.*r)+(Bn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)...
   % /Dn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3))*besselh(n,1,k1.*r)).*sin(n.*theta)+
   U=-2*pi*rp^(3)*rho1*(p1*(f1/(3*rho1^(2)*c1^(2)))-u*(f2/2));
end

