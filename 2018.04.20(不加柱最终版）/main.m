clc
clear
%设置参数
n_max=3;
scale=1e-3;
sigma=0.14;
E=7.2e10;
scale_nm=1e-9;
%毫米
a=1.5*scale;
b=1.2*scale;
rp=250*scale_nm;
%千克/立方米
rho1=1000;
rho2=2200;
rho3=1000;
rhop=1050;
%米/秒
c1=1490;
cl2=sqrt((E*(1-sigma))/(rho2*(1+sigma)*(1-2*sigma)));
ct2=sqrt(E/(2*rho2*(1+sigma)));
c3=1490;
cp=1150;
%赫兹
k_min=0;
k_max=50/scale;
omega_min=c1*k_min;
omega_max=c1*k_max;
mode=1;%解的个数 
precision=1e-5;%精确度
%**************************************************************************
%求根
omega=NaN*zeros(n_max,mode);
for n=2:n_max
    fprintf('正在计算%d阶的情况\n',n)
    [omega12,num]=feval('functionTest',omega_min,omega_max,rho1,rho2,rho3,a,b,c1,cl2,ct2,c3,n,mode,scale);
    for j=1:num
        omega1=omega12(j,1);
        omega2=omega12(j,2);
        omega(n,j)=dichotomies(omega1,omega2,n,rho1,rho2,rho3,a,b,c1,cl2,ct2,c3,precision);%利用二分法求根
    end
end
%**************************************************************************
%绘制声压图
n=1;
for i=2:n_max
    n=n+1;
        p=zeros(1200,1200);
         x=-3*a+a/200:a/200:3*a;
        y=3*a:-a/200:-3*a+a/200;
        [X,Y]=meshgrid(x,y);
        [theta,r]=cart2pol(X,Y);
        pos=find(r<=b);
        p(pos)=imag(pressureIn(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3));
%         pmax=max(p(pos));
%         pmin=min(p(pos));
%         p(pos)=(p(pos)-pmin)./(pmax-pmin);
        pos= r>b & r<a;
        p(pos)=0;
        pos=find(r>=a);
        p(pos)=imag(pressureOut(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3));
%         pmax=max(p(pos));
%         pmin=min(p(pos));
%         p(pos)=(p(pos)-pmin)./(pmax-pmin);
   
   
    figure(i)
    pcolor(x,y,p)
    colormap(jet);
%     colorbar
    shading interp
    axis off
    %axis([-1.5e-3 1.5e-3 -0.9e-3 0.9e-3])
    width=800;%图像的宽度
    height=800;%图像的高度
    left=200;%距屏幕左下角的水平距离
    bottom=200;%距屏幕左下角的垂直距离
    set(gcf,'position',[left,bottom,width,height])
%    title(['总声压场(Pa) f:' num2str(omega(i)/(2*pi))])
 index=n;
   if n>1
            save(int2str(index),'p');
   end
end


n=1;
for i=2:n_max
    n=n+1;
    
        u=zeros(600,600);
         x=-5*a+a/60:a/60:5*a;
        y=5*a:-a/60:-5*a+a/60;
        [X,Y]=meshgrid(x,y);
        [theta,r]=cart2pol(X,Y);
        pos=find(r<=b);
        u(pos)=uin(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp);
%    pmax=max(p(pos));
%  pmin=min(p(pos));
%  p(pos)=(p(pos)-pmin)./(pmax-pmin);
        pos= find(r>b & r<a);
% u(pos)=real(uke(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3));
       u(pos)=NaN;
        pos=find(r>=a);
        u(pos)=uout(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp);
%  pmax=max(p(pos));
%   pmin=min(p(pos));
% p(pos)=(p(pos)-pmin)./(pmax-pmin);
   
   
    figure(i)
    pcolor(x,y,abs(u))
    colormap(jet);
    colorbar
    shading interp
    axis off
    %axis([-1.5e-3 1.5e-3 -0.9e-3 0.9e-3])
    width=800;%图像的宽度
    height=700;%图像的高度
    left=200;%距屏幕左下角的水平距离
    bottom=200;%距屏幕左下角的垂直距离
    set(gcf,'position',[left,bottom,width,height])
   title(['总声压场(Pa) f:' num2str(omega(i)/(2*pi))])
   hold on 
   index=n;
   if n>1
            save(int2str(index),'u');
   end
   if mod(i,2)==0;
vv=zeros(250,250);   uu=vv;
         x=-5*a+a/25:a/25:5*a;
        y=5*a:-a/25:-5*a+a/25;
         [X,Y]=meshgrid(x,y);
        [theta,r]=cart2pol(X,Y);
         pos=find(r<=b);
        uu(pos)=real(uuin(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
        vv(pos)=real(vvin(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
         pos= find(r>b & r<a);
       uu(pos)=NaN;
        vv(pos)=NaN;
        pos=find(r>=a);
        uu(pos)=real(uuout(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
        vv(pos)=real(vvout(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
        U=sqrt(uu.^2+vv.^2);
        quiver(X,Y,uu./U,vv./U,0.5,'w')
    else
%          u=zeros(600,600);
%          x=-5*a+a/60:a/60:5*a;
%         y=5*a:-a/60:-5*a+a/60;
%         [X,Y]=meshgrid(x,y);
%         [theta,r]=cart2pol(X,Y);
%         pos=find(r<=b);
%         u(pos)=real(uin(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
% %    pmax=max(p(pos));
% %  pmin=min(p(pos));
% %  p(pos)=(p(pos)-pmin)./(pmax-pmin);
%         pos= find(r>b & r<a);
% % u(pos)=real(uke(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3));
%        u(pos)=NaN;
%         pos=find(r>=a);
%         u(pos)=real(uout(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
% %  pmax=max(p(pos));
% %   pmin=min(p(pos));
% % p(pos)=(p(pos)-pmin)./(pmax-pmin);
%    
%    
%     figure(i)
%     pcolor(x,y,abs(u))
%     colormap(jet);
%     colorbar
%     shading interp
%     axis off
%     %axis([-1.5e-3 1.5e-3 -0.9e-3 0.9e-3])
%     width=800;%图像的宽度
%     height=700;%图像的高度
%     left=200;%距屏幕左下角的水平距离
%     bottom=200;%距屏幕左下角的垂直距离
%     set(gcf,'position',[left,bottom,width,height])
%    title(['总声压场(Pa) f:' num2str(omega(i)/(2*pi))])
%    hold on 
%   if n>1
%             save(int2str(index),'u');
%    end
vv=zeros(250,250);   uu=vv;
         x=-5*a+a/25:a/25:5*a;
        y=5*a:-a/25:-5*a+a/25;
         [X,Y]=meshgrid(x,y);
        [theta,r]=cart2pol(X,Y);
         pos=find(r<=b);
        uu(pos)=real(uuin(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
        vv(pos)=real(vvin(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
         pos= find(r>b & r<a);
       uu(pos)=NaN;
        vv(pos)=NaN;
        pos=find(r>=a);
        uu(pos)=real(uuout(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
        vv(pos)=real(vvout(r(pos),theta(pos),n,rho1,rho2,rho3,a,b,omega(i),c1,cl2,ct2,c3,rp,rhop,cp));
        U=sqrt(uu.^2+vv.^2);
        quiver(X,Y,uu./U,vv./U,0.5,'w')
    end
end

