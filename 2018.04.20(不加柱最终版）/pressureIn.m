function p3=pressureIn(r,theta,n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)
i=sqrt(-1);
k3=omega/c3;
if n==0
    epsilon=1;
else
    epsilon=2;
end
p3=i^(n)*epsilon*qn(n,rho1,rho2,rho3,a,b,omega,c1,cl2,ct2,c3)*besselj(n,k3.*r).*cos(n.*theta);
end

