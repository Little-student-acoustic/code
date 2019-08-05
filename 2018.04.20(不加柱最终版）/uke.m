function u=uke(r,theta,n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3,rp,rhop,cp)
i=sqrt(-1);
syms x1  x2
dbesselj=matlabFunction(diff(besselj(n,x1),x1));
dbessely=matlabFunction(diff(bessely(n,x2),x2));
k_l2=omega/cl2;
k_t2=omega/ct2;
if n==0
    epsilon=1;
else
    epsilon=2;
end
u=i^(n)*epsilon*((-n/r)*(gn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*besselj(n,k_l2.*r)+hn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*bessely(n,k_l2.*r))-k_t2*(ln(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbesselj(k_t2.*r)+mn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbessely(k_t2.*r))).*sin(n.*theta)+i^(n)*epsilon*cos(n.*theta)*(k_l2*gn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbesselj(k_t2.*r)+hn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*dbessely(k_t2.*r)+(n/r)*(ln(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*besselj(n,k_l2.*r)+mn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*bessely(n,k_l2.*r)));
end

