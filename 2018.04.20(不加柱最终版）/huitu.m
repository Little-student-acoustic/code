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
    width=500;%ͼ��Ŀ��
    height=400;%ͼ��ĸ߶�
    left=200;%����Ļ���½ǵ�ˮƽ����
    bottom=200;%����Ļ���½ǵĴ�ֱ����
    set(gcf,'position',[left,bottom,width,height])
   % title('����ѹ��(Pa)')
%}
 save(int2str(i),'x','y','p','part');
end
