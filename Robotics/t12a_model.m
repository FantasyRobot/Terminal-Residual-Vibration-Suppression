classdef t12a_model < handle
    properties
        a1;
        a3;
        d0;
        d1;
        d2;
        d3;
        d4; 
        d5;
    end
    methods
        function obj = t12a_model(para_name)
            switch para_name
                case 'test'
             %% for calib data
                    obj.a1=0.1;
                    obj.a3=0.135;
                    obj.d0=0;
                    obj.d1=0.615;
                    obj.d2=0.705;
                    obj.d3=0;
                    obj.d4=0.3775*2;                     
                    obj.d5=0.085;  
                otherwise
            %% t12a normal param
                    obj.a1=0.185;
                    obj.a3=0.2;
                    obj.d0=0;
                    obj.d1=0;
                    obj.d2=0.6;
                    obj.d3=0;
                    obj.d4=0.65; 
                    obj.d5=0.105; 
            end         
        end
       
        %% 论文坐标系定义DH标准模型正解
        function [ R, P, C ] = kine_normal( obj, theta)
            R01=[cos(theta(1,1)), -sin(theta(1,1)),0;...
                 sin(theta(1,1)),cos(theta(1,1)),0;...
                 0,0,1];
            R12=[cos(theta(2,1)), 0, sin(theta(2,1));...
                 0,1,0;...
                 -sin(theta(2,1)),0,cos(theta(2,1))];
            R23=[cos(theta(3,1)), 0, sin(theta(3,1));...
                 0,1,0;...
                 -sin(theta(3,1)),0,cos(theta(3,1))];
            R34=[1, 0, 0;...
                 0,cos(theta(4,1)),-sin(theta(4,1));...
                 0,sin(theta(4,1)),cos(theta(4,1))];
            R45=[cos(theta(5,1)), 0, sin(theta(5,1));...
                 0,1,0;...
                 -sin(theta(5,1)),0,cos(theta(5,1))];
            R56=[1, 0, 0;...
                 0,cos(theta(6,1)),-sin(theta(6,1));...
                 0,sin(theta(6,1)),cos(theta(6,1))];
            R=R01*R12*R23*R34*R45*R56;
%             Pw=[0;0;obj.d0+obj.d1]+R01*[obj.a1;0;0]+(R01*R12)*[0;0;obj.d2]+(R01*R12*R23)*[obj.d3;0;obj.a3]+(R01*R12*R23*R34)*[obj.d4;0;0]+(R01*R12*R23*R34*R45)*[obj.d5;0;0];
            Pw=[0;0;obj.d0+obj.d1]+R01*[obj.a1;0;0]+(R01*R12)*[0;0;obj.d2]+(R01*R12*R23)*[obj.d3+obj.d4;0;obj.a3];

            P=Pw+R*[obj.d5;0;0];
            C =0;
        end
        %% 论文坐标系定义DH标准模型反解
        function [Q]= ikine_normal( obj, P, R)
            pwx=P(1,1);pwy=P(2,1);pwz=P(3,1);
            % 求解关节1
            Q(1,1) = atan(pwy/pwx);

            % 求解关节2
            aa = cos(Q(1,1))*pwx+sin(Q(1,1))*pwy-obj.a1;
            a = -2*aa*obj.d2;
            bb = pwz-(obj.d0+obj.d1);
            b = -2*bb*obj.d2;
            c = aa^2+bb^2+obj.d2^2-obj.a3^2-(obj.d3+obj.d4)^2;     
            a^2+b^2-c^2
            Q(2,1) = 2*atan((-2*a+(a^2+b^2-c^2)^0.5)/(-b+c));
            %Q(2,1) = 2*atan((-2*a-(a^2+b^2-c^2)^0.5)/(-b+c));
                       
        end
        %% robotics 库DH标准模型正解
        function [robot] = robot_tool(obj)
            theta1 = 0;     D1 = obj.d1; A1 = obj.a1; alpha1 = -pi/2;  offset1 = 0;
            theta2 = -pi/2; D2 = 0;     A2 = obj.d2;   alpha2 = 0;      offset2 = pi/2;
            theta3 = 0;     D3 = 0;     A3 = obj.a3;   alpha3 = -pi/2;  offset3 = 0;
            theta4 = 0;     D4 = obj.d4;  A4 = 0;     alpha4 = pi/2;   offset4 = 0;
            theta5 = 0;      D5 = 0;    A5 = 0;     alpha5 = -pi/2;  offset5 = 0;
            theta6 = 0;     D6 =obj.d5; A6 = 0;     alpha6 = 0;      offset6 = 0;
            L(1) = Link([theta1, D1, A1, alpha1, offset1], 'standard');
            L(2) = Link([theta2, D2, A2, alpha2, offset2], 'standard');
            L(3) = Link([theta3, D3, A3, alpha3, offset3], 'standard');
            L(4) = Link([theta4, D4, A4, alpha4, offset4], 'standard');
            L(5) = Link([theta5, D5, A5, alpha5, offset5], 'standard');
            L(6) = Link([theta6, D6, A6, alpha6, offset6], 'standard');
            %t12a_mdh = SerialLink([R1,R2,R3,R4,R5,R6]);
            robot = SerialLink(L);      
        end

        %% 论文坐标系定义DH关节结合误差模型正解
        function [ R, P] = kine_abs(obj, theta, abs_err)
            R01=rotx(abs_err(1,1))*roty(abs_err(1,2))*...
                [cos(theta(1,1)-abs_err(1,3)), -sin(theta(1,1)-abs_err(1,3)),0;...
                 sin(theta(1,1)-abs_err(1,3)),cos(theta(1,1)-abs_err(1,3)),0;...
                 0,0,1];
            R12=rotx(abs_err(2,1))*rotz(abs_err(2,2))*...
                [cos(theta(2,1)-abs_err(2,3)), 0, sin(theta(2,1)-abs_err(2,3));...
                 0,1,0;...
                 -sin(theta(2,1)-abs_err(2,3)),0,cos(theta(2,1)-abs_err(2,3))];
            R23=rotx(abs_err(3,1))*rotz(abs_err(3,2))*...
                [cos(theta(3,1)-abs_err(3,3)), 0, sin(theta(3,1)-abs_err(3,3));...
                 0,1,0;...
                 -sin(theta(3,1)-abs_err(3,3)),0,cos(theta(3,1)-abs_err(3,3))];
            R34=roty(abs_err(4,1))*rotz(abs_err(4,2))*...
                [1, 0, 0;...
                 0,cos(theta(4,1)-abs_err(4,3)),-sin(theta(4,1)-abs_err(4,3));...
                 0,sin(theta(4,1)-abs_err(4,3)),cos(theta(4,1)-abs_err(4,3))];
            R45=rotx(abs_err(5,1))*rotz(abs_err(5,2))*...
                [cos(theta(5,1)-abs_err(5,3)), 0, sin(theta(5,1)-abs_err(5,3));...
                 0,1,0;...
                 -sin(theta(5,1)-abs_err(5,3)),0,cos(theta(5,1)-abs_err(5,3))] 
            R56=roty(abs_err(6,1))*rotz(abs_err(6,2))*...
                [1, 0, 0;...
                 0,cos(theta(6,1)-abs_err(6,3)),-sin(theta(6,1)-abs_err(6,3));...
                 0,sin(theta(6,1)-abs_err(6,3)),cos(theta(6,1)-abs_err(6,3))]; 
            R=R01*R12*R23*R34*R45*R56;
            P=[abs_err(1,4);abs_err(1,5);abs_err(1,6)+obj.d0+obj.d1]+...
               R01*[abs_err(2,4)+obj.a1;abs_err(2,5);abs_err(2,6)]+...
               R01*R12*[abs_err(3,4);abs_err(3,5);abs_err(3,6)+obj.d2]+...
               R01*R12*R23*[abs_err(4,4)+obj.d3;abs_err(4,5);abs_err(4,6)+obj.a3]+...
               R01*R12*R23*R34*[abs_err(4,4)+obj.d4;abs_err(4,5);abs_err(4,6)]+...
               R01*R12*R23*R34*R45*[abs_err(5,4)+obj.d5;abs_err(5,5);abs_err(5,6)]+...
               R*[0;0;0];
        end
        
    end
end