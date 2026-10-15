function [P,R]=T12_fkine(theta)
%% 杆长参数
obj.a1=0.185; % L1
obj.d2=0.6;   % L2
obj.a3=0.2;   % L3
obj.d4=0.65;  % L4
obj.d5=0.105; % L5
obj.d0=0;
obj.d1=0;
obj.d3=0;

%% 关节角度
% theta = [-pi/4;-pi/4;-pi/4;-pi/4;-pi/4;-pi/4];
% theta = [pi/4;pi/4;pi/4;pi/4;pi/4;pi/4];
%  theta = [0;0;0;0;0;0];
theta = theta.*pi/180;

R01=[cos(theta(1,1)), -sin(theta(1,1)),0;...
     sin(theta(1,1)),cos(theta(1,1)),0;...
     0,0,1];
R12=[cos(-theta(2,1)), 0, sin(-theta(2,1));...
     0,1,0;...
     -sin(-theta(2,1)),0,cos(-theta(2,1))];
R23=[cos(-theta(3,1)), 0, sin(-theta(3,1));...
     0,1,0;...
     -sin(-theta(3,1)),0,cos(-theta(3,1))];
R34=[1, 0, 0;...
     0,cos(theta(4,1)),-sin(theta(4,1));...
     0,sin(theta(4,1)),cos(theta(4,1))];
R45=[cos(-theta(5,1)), 0, sin(-theta(5,1));...
     0,1,0;...
     -sin(-theta(5,1)),0,cos(-theta(5,1))];
R56=[1, 0, 0;...
     0,cos(theta(6,1)),-sin(theta(6,1));...
     0,sin(theta(6,1)),cos(theta(6,1))];
%% 关节3位置P3
P3=[0;0;obj.d0+obj.d1]+R01*[obj.a1;0;0]+(R01*R12)*[0;0;obj.d2]; 

%% X3与Z4交点处位置P3_1
P3_1=[0;0;obj.d0+obj.d1]+R01*[obj.a1;0;0]+(R01*R12)*[0;0;obj.d2]+(R01*R12*R23)*[obj.d3;0;obj.a3];

%% 末端位置P
R=R01*R12*R23*R34*R45*R56;
Pw=[0;0;obj.d0+obj.d1]+R01*[obj.a1;0;0]+(R01*R12)*[0;0;obj.d2]+(R01*R12*R23)*[obj.d3+obj.d4;0;obj.a3];
P=Pw+R*[obj.d5;0;0];
end
%% DH法建模
% theta1 = 0;     D1 = obj.d1; A1 = obj.a1; alpha1 = -pi/2;  offset1 = 0;
% theta2 = -pi/2; D2 = 0;     A2 = obj.d2;   alpha2 = 0;      offset2 = pi/2;
% theta3 = 0;     D3 = 0;     A3 = obj.a3;   alpha3 = -pi/2;  offset3 = 0;
% theta4 = 0;     D4 = obj.d4;  A4 = 0;     alpha4 = pi/2;   offset4 = 0;
% theta5 = 0;      D5 = 0;    A5 = 0;     alpha5 = -pi/2;  offset5 = 0;
% theta6 = 0;     D6 =obj.d5; A6 = 0;     alpha6 = 0;      offset6 = 0;
% L(1) = Link([theta1, D1, A1, alpha1, offset1], 'standard');
% L(2) = Link([theta2, D2, A2, alpha2, offset2], 'standard');
% L(3) = Link([theta3, D3, A3, alpha3, offset3], 'standard');
% L(4) = Link([theta4, D4, A4, alpha4, offset4], 'standard');
% L(5) = Link([theta5, D5, A5, alpha5, offset5], 'standard');
% L(6) = Link([theta6, D6, A6, alpha6, offset6], 'standard');
% %t12a_mdh = SerialLink([R1,R2,R3,R4,R5,R6]);
% robot = SerialLink(L);  