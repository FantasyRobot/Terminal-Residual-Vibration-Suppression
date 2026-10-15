clear all;clc;close all;
import casadi.*;

% [xx,ww,vv]=chebpts(12);% cheb pnts and weights
% D=ChebTest.getDifferentiationMatrix(xx,vv);% differential matrix

jnt2 = load("jnt2.txt");
jnt3 = load("jnt3.txt");
jnt5 = load("jnt5.txt");

figure (1)
subplot (2,1,1)
x=[1,12,200,453,800,910, 940,949];
y=[0,-0.004271,-0.08778,-0.07368,-0.09159,-0.01638,-0.002271,0];
xq = linspace(1, 949, 948);
vel2=zeros(948,1);
for i =1:length(jnt2)-1
    vel2(i,1) = jnt2(i+1,3)-jnt2(i,3);
end
p2 = pchip(x,y,xq);
plot(x,y,'o',xq,p2,'-',xq,vel2,'-.')
legend('Sample Points','pchip','original','Location','SouthEast');
hold on;

subplot (2,1,2)
x=[1,17,117,501,747, 900,941,949];
y=[0,0.01371,0.06385,0.08868,0.1193,0.1504,0.01144,0];
xq = linspace(1, 949, 948);
vel3=zeros(948,1);
for i =1:length(jnt3)-1
    vel3(i,1) = jnt3(i+1,3)-jnt3(i,3);
end
p3 = pchip(x,y,xq);
plot(x,y,'o',xq,p3,'-',xq,vel3,'-.')
legend('Sample Points','pchip','original','Location','SouthEast');

figure(2)
posRes2=zeros(948,1);
posRes2(1,1)=jnt2(1,3);
for i =1:length(jnt2)-1
    posRes2(i+1,1) = posRes2(i,1)+p2(1,i);
end
plot(posRes2);
hold on;
posRes3=zeros(948,1);
posRes3(1,1)=jnt3(1,3);
for i =1:length(jnt3)-1
    posRes3(i+1,1) = posRes3(i,1)+p3(1,i);
end
plot(posRes3);


fid = fopen('pathInfo.txt', 'wt');
for i = 1:length(posRes2)
    fprintf(fid,'%f,%f,%f,%f,%f,%f\n',0,posRes2(i,1),posRes3(i,1),0,jnt5(i,3),0);
end
fclose(fid);