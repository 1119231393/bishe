function [G,gabout]=Gabor_hy(I,Sx,Sy,fi,miu,sigma)
% %miu是中心频率，fi是滤波器的方向，sigma是高斯窗的标准差,Sx,Sy为空间域像素的位置，
% x=-fix(Sx):fix(Sx);%Gabor变换核函数的窗口长度,
% y=-fix(Sy):fix(Sy);
% [x,y]=meshgrid(x,y);
% xPrime=x*cos(fi)+y*sin(fi);
% yPrime=y*cos(fi)-x*sin(fi);
% %Gabor变换核函数
% %Gabor=-((1/(2*pi*sigma.^2)).*exp(-.5*(xPrime.^2+yPrime.^2)/sigma.^2).*(exp(j*f*xPrime)-exp(-(f*sigma)^2/2)));
% % Gabor=((1/(2*pi*sigma.^2)).*exp(-.5*(xPrime.^2+yPrime.^2)/sigma.^2).*exp(2*pi*sqrt(-1)*(f*xPrime*cos(fi)+f*yPrime*sin(fi))));
% Gabor=((1/(2*pi*sigma.^2)).*exp(-.5*(xPrime.^2+yPrime.^2)/sigma.^2).*cos(2*pi*sqrt(-1)*(f*xPrime*cos(fi)+f*yPrime*sin(fi))));%实部
if isa(I,'double')~=1
    I=double(I);
end
GRealSum=0;
for x = -fix(Sx):fix(Sx)
    for y = -fix(Sy):fix(Sy)
        G(fix(Sx)+x+1,fix(Sy)+y+1) = ((1/(2*pi*sigma.^2)).*exp(-.5*(x.^2+y.^2)/sigma.^2).*exp(2*pi*sqrt(-1)*(miu*x*cos(fi)+miu*y*sin(fi))));
        GReal(fix(Sx)+x+1,fix(Sy)+y+1) = ((1/(2*pi*sigma.^2)).*exp(-.5*(x.^2+y.^2)/sigma.^2).*cos(2*pi*sqrt(-1)*(miu*x*cos(fi)+miu*y*sin(fi))));
        GRealSum=GRealSum+GReal(fix(Sx)+x+1,fix(Sy)+y+1);
    end
end
for x = -fix(Sx):fix(Sx)
    for y = -fix(Sy):fix(Sy)
        GRealDeAvg(fix(Sx)+x+1,fix(Sy)+y+1)=GReal(fix(Sx)+x+1,fix(Sy)+y+1)-GRealSum/(2*Sx+1)^2;
    end
end
Regabout = conv2(I,double(GRealDeAvg),'same');
%gabout = sqrt(Imgabout.*Imgabout+Regabout.*Regabout);
gabout = Regabout;
