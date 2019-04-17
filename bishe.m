close all; 
clear; 
clc; 
PathRoot='C:\Users\hasee\Desktop\graduation project\palmbase';
list=dir(PathRoot);
fileNum=size(list); 
%for k=3:fileNum
for k=3:fileNum
	disp(list(k).name)  % 这就是文件名，如果有子文件夹，则也包含在里面。
%f1=imread('1479305870650c.jpg');
f1=imread(list(k).name);
f1=imresize(f1,0.33);
[m,n]=size(f1);
n=n/3;
figure(1);imshow(f1);
hold on
xA = 2*n/5;
yA = m/3;
xB = 3*n/5;
yB = yA - (xB - xA)*tand(30);
xC = 2*n/5;
yC = 4*m/5;
xD = n;
yD = yC - (xD - xC)*tand(30);
plot([xA,xB],[yA,yB],'Color','b','LineWidth',5);%AB连线
plot([xC,xD],[yC,yD],'Color','b','LineWidth',5);%CD连线
plot(xB,yB,'r.','markersize',20);%画红点
hold off;
%%旋转
f1=imrotate(f1,-30,'bilinear','crop');
figure(2);imshow(f1);
hold on
x2=sqrt((xB-n/2)*(xB-n/2)+(yB-m/2)*(yB-m/2))*cosd(atand(abs((yB-m/2)/(xB-n/2)))-30); 
y2=sqrt((xB-n/2)*(xB-n/2)+(yB-m/2)*(yB-m/2))*sind(atand(abs((yB-m/2)/(xB-n/2)))-30); 
xB2=n/2+x2;
yB2=m/2-y2;
plot(xB2,yB2,'r.','markersize',20);%画旋转后B点



x2=sqrt((xC-n/2)*(xC-n/2)+(yC-m/2)*(yC-m/2))*cosd(atand((m/2-yC)/(xC-n/2))-30); 
y2=sqrt((xC-n/2)*(xC-n/2)+(yC-m/2)*(yC-m/2))*sind(atand((m/2-yC)/(xC-n/2))-30); 
xC2=n/2-x2;
yC2=m/2+y2;
plot(xC2,yC2,'b.','markersize',20);%画旋转后C点


l=yC2-yB2;
L= 13*l/20;
xE=xB2-3*l/4;
yE=yB2+3*l/20;
plot(xE,yE,'g.','markersize',20);
%画旋转后ROI左上角，l 为手掌上下边界之间的间距；LROI 为 ROI 边长；E 为 ROI 左上角的顶点。
d=l/10;%d 为偏移量，控制区域范围调整，d=l/10。候选区域宽设为 d，高设为 3d。
xF = xB2 - 2*d;
yF = yB2 - d;
xG = xB2 - 4*d;
yG = yC2 - 2*d;
% for i=yC2-l/10:yC2+l/10 %G为左上角的上候选区，候选区域方向编号
%    for j=xB2-3*l/5:xB2-l/5
yg2=yC2-l/10;
xg2=xB2-3*l/5;
plot(xg2,yg2,'g.','markersize',20);
plot(xF,yF,'w.','markersize',20);
plot(xG,yG,'w.','markersize',20);
hold off;
theta1=0*pi/3;
theta2=1*pi/3;
theta3=2*pi/3;
% f=0.88;
% sigma=2.6;

f=0.0916;
sigma=5.6179;
Sx=5;%窗口长度
Sy=5;
% Gabor1=real(Gabor_hy(Sx,Sy,f,theta1,sigma));
% Gabor2=real(Gabor_hy(Sx,Sy,f,theta2,sigma));
% Gabor3=real(Gabor_hy(Sx,Sy,f,theta3,sigma));
Gabor1=Gabor_hy(Sx,Sy,f,theta1,sigma);
Gabor2=Gabor_hy(Sx,Sy,f,theta2,sigma);
Gabor3=Gabor_hy(Sx,Sy,f,theta3,sigma);
%Gabor1junzhi=Gabor1-(Gabor1+Gabor2+Gabor3)/3;
%Gabor2junzhi=Gabor2-(Gabor1+Gabor2+Gabor3)/3;
%Gabor3junzhi=Gabor3-(Gabor1+Gabor2+Gabor3)/3;

Gabor1junzhi=Gabor1-sum(Gabor1(:))/(2*Sx+1)^2;
Gabor2junzhi=Gabor2-sum(Gabor2(:))/(2*Sx+1)^2;
Gabor3junzhi=Gabor3-sum(Gabor3(:))/(2*Sx+1)^2;
hold on
plot(1,1,'r.','markersize',20);%画红点
hold off;
Regabout1=conv2(double(rgb2gray(f1)),double(Gabor1junzhi),'same');
Regabout2=conv2(double(rgb2gray(f1)),double(Gabor2junzhi),'same');
Regabout3=conv2(double(rgb2gray(f1)),double(Gabor3junzhi),'same');
%三个滤波结果的最小值的绝对值作为该点的能量，并记录最小值对应的方向编号，以上候选区域为例：
Regabout=zeros(m,n);
Regabout0=Regabout;
direction=zeros(m,n);

% for i=182:182+3*d+1 %F为左上角的上候选区，候选区域方向编号
%    for j=231:231+d+1
for i=floor(yB2-l/10):yB2+l/10 %F为左上角的上候选区，候选区域能量
    for j=floor(xB2-l/3):xB2

        if(Regabout1(i,j)<=Regabout2(i,j))
            min=Regabout1(i,j);
            minDirection=1;
        else
            min=Regabout2(i,j);
            minDirection=2;
        end
        if(min>Regabout3(i,j))
            min=Regabout3(i,j);
            minDirection=3;
        end
        fprintf('min=%f,minDirection=%f;',min, minDirection);
        Regabout(i,j)=abs(min);
        direction(i,j)=minDirection;
    end
end   

% for i=182:182+3*d+1 %F为左上角的上候选区，候选区域方向编号
%    for j=231:231+d+1
for i=floor(yB2-l/10):yB2+l/10 %F为左上角的上候选区，候选区域能量
    for j=floor(xB2-l/3):xB2
        if(Regabout1(i,j)<Regabout2(i,j))
            min=Regabout1(i,j);
            minDirection=1;
        else
            min=Regabout2(i,j);
            minDirection=2;
        end
        if(min>Regabout3(i,j))
            min=Regabout3(i,j);
            minDirection=3;
        end
        fprintf('min=%f,minDirection=%f;',min, minDirection);
        if(minDirection==1)
        Regabout0(i,j)=abs(min);
        else
        Regabout0(i,j)=0;
        end
    end
end      

% for i=412:412+3*d+1 %G为左上角的上候选区，候选区域方向编号
%    for j=180:180+d+1
for i=floor(yC2-l/10):yC2+l/10 %G为左上角的上候选区，候选区域方向编号
   for j=floor(xB2-3*l/5):xB2-l/5
        if(Regabout1(i,j)<Regabout2(i,j))
            min=Regabout1(i,j);
            minDirection=1;
        else
            min=Regabout2(i,j);
            minDirection=2;
        end
        if(min>Regabout3(i,j))
            min=Regabout3(i,j);
            minDirection=3;
        end
        fprintf('min=%f,minDirection=%f;',min, minDirection);
        Regabout(i,j)=abs(min);
        direction(i,j)=minDirection;
    end
end   

% for i=412:412+3*d+1 %G为左上角的上候选区，候选区域方向编号
%    for j=180:180+d+1
for i=floor(yC2-l/10):yC2+l/10 %G为左上角的上候选区，候选区域方向编号
   for j=floor(xB2-3*l/5):xB2-l/5
        if(Regabout1(i,j)<Regabout2(i,j))
            min=Regabout1(i,j);
            minDirection=1;
        else
            min=Regabout2(i,j);
            minDirection=2;
        end
        if(min>Regabout3(i,j))
            min=Regabout3(i,j);
            minDirection=3;
        end
        fprintf('min=%f,minDirection=%f;',min, minDirection);
        if(minDirection==1)
        Regabout0(i,j)=abs(min);
        else
        Regabout0(i,j)=0;
        end
    end
end     

set(0,'defaultFigurePosition',[100,100,1200,450]);
set(0,'defaultFigureColor',[1 1 1]);
figure(3),
subplot(141),imshow(f1);
subplot(142),imshow(Regabout1);
subplot(143),imshow(Regabout2);
subplot(144),imshow(Regabout3);
figure(4),
imshow(Regabout);title('选区域能量');
figure(5),
imshow(Regabout0);
% 求水平投影
for x=1:m
    S(x)=sum(Regabout0(x,:));
end
x=1:m;
figure(6),
plot(x,S(x));
title('水平投影');
%
%[a,b]=findpeaks(S(x));%a,b 每个峰值的数值及其横坐标
mu=max(S(182<x & x<182+3*d+1));
yu=182;
for u=182:182+3*d+1
    if(S(u)>mu/2)
        yu=u;
    end
end
md=max(S(412<x & x<412+3*d+1));%其中0.00001是精度。412:412+3*d+1
yd=463;
for u=412:412+3*d+1
    if(S(u)>md/2)
        yd=u;
        break;
    end
end
figure(7),
imshow(f1);
hold on
%
plot([0,n],[yB2,yB2],'Color','b','LineWidth',2);
plot([0,n],[yC2,yC2],'Color','b','LineWidth',2);
%水平积分投影后实际边界，上边界线的 y 轴坐标 yu 为 pu 中的最大值，即 yu=max(pu)；类似的，下边界线的 y 轴 坐 标 yd 为 pd 中 的 最 小 值 ， 即yd=min(pd)。
plot([0,n],[yu,yu],'Color','r','LineWidth',2);
plot([0,n],[yd,yd],'Color','r','LineWidth',2);
%根据实际边界找ROI
l2=yd-yu;
L2=7*l2/10;
xE2=xB2-3*l2/4;
yE2=yu+l2/8;
plot(xE2,yE2,'r.','markersize',20);%画红点
saveas(gcf,['C:\Users\hasee\Desktop\graduation project\preprocessresult\',list(k).name,'.jpg']);
hold off;
% 
% %%
% %截取
% RGB=f1;
% RGB1=imcrop(RGB,[xG,yG,d,3*d]);
% figure(22),
% imshow(RGB1);
% saveas(RGB1,['C:\Users\hasee\Desktop\graduation project\preprocessresult\',RGB1,'.jpg']);
% %imwrite(RGB1,['C:\Users\hasee\Desktop\graduation project\preprocessresult\.RGB1jpg']);
% 


 end
