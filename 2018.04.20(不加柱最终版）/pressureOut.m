function p1=pressureOut(r,theta,n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)
i=sqrt(-1);
k1=omega/c1;
if n==0
    epsilon=1;
else
    epsilon=2;
end
p1=i^(n)*epsilon*(besselj(n,k1.*r)+(Bn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)...
    /Dn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3))*besselh(n,1,k1.*r)).*cos(n.*theta);
end

