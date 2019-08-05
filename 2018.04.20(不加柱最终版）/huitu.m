clc
clear
for i=2:7

   load(int2str(i));

part=abs(p(650,:));
subplot(2,3,i-1);
plot(x,part,'b');
%hold on

   %{
figure(i)
pcolor(x,y,p)
colormap(jet);
   % colorbar
  shading interp
   axis off
    width=500;%图像的宽度
    height=400;%图像的高度
    left=200;%距屏幕左下角的水平距离
    bottom=200;%距屏幕左下角的垂直距离
    set(gcf,'position',[left,bottom,width,height])
   % title('总声压场(Pa)')
%}
 save(int2str(i),'x','y','p','part');
end
