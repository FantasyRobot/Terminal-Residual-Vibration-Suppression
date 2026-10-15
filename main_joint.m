clear all;clc;close all;

addpath(genpath('./PseudoOpt'));
addpath(genpath('./Robotics'));
addpath(genpath('./Motion'));
addpath(genpath('./Data/initData'));

jnt2 = load("jnt2.txt");
jnt3 = load("jnt3.txt");
jnt5 = load("jnt5.txt");
tcp = load("tcp.txt")
ps=PseudoOptimal;
ps.npts=16;
ps.nS=12;
ps.nU=6;
% maxiter=500;%过多的优化次数使轨迹变化过大
maxiter=1;
dispLev=5;

prob.rob = RobViz('M20iA');

[prob.c,prob.r,prob.obs,...
    prob.selfmap,prob.wallmap,prob.obsmap] =...
    boundBall(prob);

% prob.x_init=zeros(6,ps.npts);

ps=ps.timeOptimalInit(@(Xc,Uc,D,nS,nU,scale,P)conEq(Xc,Uc,D,nS,nU,scale,P,prob),...
                      @(Xc,Uc,tc,D,nS,nU,ww,scale)costFcn_topt(Xc,Uc,tc,D,nS,nU,ww,scale,prob),...
                      @(Xc,Uc,D,scale)conIneq(Xc,Uc,D,scale,prob),...
                      maxiter,...
                      dispLev);                

[prob.lb,prob.ub,prob.bnds]=bounds(ps.npts,zeros(1,6),zeros(1,6),prob);
velbnd=prob.bnds(:,3);
trqbnd=prob.bnds(:,4);
dtrqbnd=prob.bnds(:,5);

state_record=[];control_record=[];time_record=[];

prob.init = [0, -25*pi/180, 54*pi/180, 0, 60*pi/180, 0];% joint space initial position
prob.target = [0,50*pi/180, -32*pi/180, 0, 72*pi/180, 0];% joint space target position
prob.targetPos = prob.rob.tcpPos(prob.target);
disp(['Initial pos [',num2str(prob.init(:).'),']']);
disp(['Target pos [',num2str(prob.targetPos(:).'),']']);

%prob.tf=20*max(abs(prob.target(:)-prob.init(:))./velbnd(:));
prob.dt = 0.004;
prob.tf=length(jnt2)*prob.dt;
time = 0:prob.dt:prob.tf;
prob.time =time(:,1:end-1);

%% generate initial motion
ps.sGuess=zeros(length(prob.time),ps.nS);
sAcc=zeros(length(prob.time),ps.nS/2);
prob.u_init=zeros(length(prob.time),ps.nU);

data_length = length(prob.time);
jnt_init(:,1) = zeros(data_length,1);
jnt_init(:,2) = -jnt2(1:data_length,3);
jnt_init(:,3) = -jnt3(1:data_length,3);
jnt_init(:,4) = zeros(data_length,1);
jnt_init(:,5) = -jnt5(1:data_length,3);
jnt_init(:,6) = zeros(data_length,1);
vel_init=zeros(data_length,6);
acc_init=zeros(data_length,6);
jerk_init=zeros(data_length,6);
for j=1:6
    for i=1:data_length-1
        vel_init(i,j) = (jnt_init(i+1,j)-jnt_init(i,j));        
    end
end
for i=1:data_length-2
    for j=1:6
        acc_init(i,j) = vel_init(i+1,j)-vel_init(i,j);
    end
end
for i=1:data_length-3
    for j=1:6
        jerk_init(i,j) = acc_init(i+1,j)-acc_init(i,j);
    end
end
%% 方法1：生成初始路径
% for i=1:ps.nS/2
%     [pos,vel,acc] = motiongen(prob.init(i),prob.target(i),prob.dt,prob.tf);
%     ps.sGuess(:,i) = pos(1:end-1);
%     ps.sGuess(:,6+i) = vel(1:end-1);
%     sAcc(:,i) = acc(1:end-1);
% end
% 
% for j=1:size(sAcc,1)
%      prob.u_init(j,:)= prob.rob.rne_s(ps.sGuess(j,1:6),ps.sGuess(j,7:12),sAcc(j,:));
% end
    
%% 方法1 end

%% 方法2：加载初始路径

ps.sGuess(:,1) = zeros(data_length,1);
ps.sGuess(:,2) = -jnt2(1:data_length,3);
ps.sGuess(:,3) = -jnt3(1:data_length,3);
ps.sGuess(:,4) = zeros(data_length,1);
ps.sGuess(:,5) = -jnt5(1:data_length,3);
ps.sGuess(:,6) = zeros(data_length,1);


for j=1:6
    for i=1:data_length-1     
        ps.sGuess(i,6+j) = vel_init(i,j);
    end
end
for i=1:data_length-2
    for j=1:6
        sAcc(i,j) = acc_init(i,j);
    end
end
%% 控制量初值，获得关节扭矩
for j=1:size(sAcc,1)
     prob.u_init(j,:)=prob.rob.rne_s(ps.sGuess(j,1:6),ps.sGuess(j,7:12),sAcc(j,:));
end
%% 方法2 end

[prob.lb,prob.ub,prob.bnds]=bounds(ps.npts,prob.init,prob.targetPos,prob);

%% solve optimal control
[opt.Xc,opt.Uc,opt.Tc,opt.Xopt,opt.Uopt,opt.Topt] = ...
    ps.timeOptimalOptimize(prob.lb,prob.ub,prob.time,prob.u_init);

%% plot recorded results
state_record = [opt.Xopt];
control_record = [opt.Uopt];
if isempty(time_record)
    time_record = [time_record;opt.Topt(:)];
else
    time_record = [time_record;opt.Topt(:)+time_record(end)+prob.dt];
end

% plot result
ratioVel=zeros(length(state_record(:,1)),6);
ratioPos=zeros(length(state_record(:,1)),6);
ratioAcc=zeros(length(state_record(:,1)),6);
ratioTrq=zeros(length(control_record(:,1)),6);
ratioJerk=zeros(length(control_record(:,1)),6);
dtrq=zeros(length(control_record(:,1)),6);
ratioDtrq=zeros(length(dtrq(:,1)),6);

for j=1:6
    ratioPos(:,j)=state_record(:,j);
    ratioVel(:,j)=state_record(:,6+j);
    % ratioTrq(:,j)=control_record(:,j)/trqbnd(j);
 
    for i=1:length(ratioVel(:,j))-1
    	ratioAcc(i,j)= ratioVel(i+1,j)-ratioVel(i,j);
    end
    for i=1:length(ratioAcc(:,j))-1
    	ratioJerk(i,j)= ratioAcc(i+1,j)-ratioAcc(i,j);
    end
    ratioTrq(:,j)=control_record(:,j);
    dtrq(:,j)=gradient(control_record(:,j))./gradient(time_record(:));
    ratioDtrq(:,j)=dtrq(:,j)/dtrqbnd(j);
end

figure(1);
clf;
subplot(2,2,1);
pos=zeros(length(state_record(:,1)),6);
for j=1:6
    pos(1,j)=jnt_init(1,j);
end
for j=1:6
    ratioVel(:,j)=state_record(:,6+j);
%     pos(1,j)=ratioPos(1,j);
end

for j = 1:6
    for i=1:length(ratioVel(:,j))-1
        pos(i+1,j)=pos(i,j)+ratioVel(i,j);
    end
end

for j=1:6
     plot(time_record,pos(:,j),'linewidth',2);
%     plot(time_record,ratioPos(:,j).*(180/pi),'linewidth',2);
    hold on;
end
legend('关节1','关节2','关节3','关节4','关节5','关节6');
grid on;xlim([time_record(1),time_record(end)]);
title('(a) 优化后位置曲线');
xlabel('时间 T(s)');ylabel('角度 \theta(rad)');

subplot(2,2,2);
for j=1:6
    C(j)=plot(time_record,ratioVel(:,j),'linewidth',2);
    hold on;
end
grid on;xlim([time_record(1),time_record(end)]);
title('(b) 优化后速度曲线');
xlabel('时间 T(s)');ylabel('速度 v(rad/s)');

subplot(2,2,3);
for j=1:6
    D(j)=plot(time_record,ratioAcc(:,j),'linewidth',2);
    hold on;
end
grid on;xlim([time_record(1),time_record(end)]);
title('(c) 优化后加速度率曲线');
xlabel('时间 T(s)');ylabel('加速度 a(rad/s^2) ');

subplot(2,2,4);
for j=1:6
    plot(time_record,ratioJerk(:,j),'linewidth',2);
    hold on;
end
grid on;xlim([time_record(1),time_record(end)]);
title('(d) 优化后Jerk率曲线');
xlabel('时间 T(s)');ylabel('jerk j(rad/s^3) ');

figure (2)
%位置曲线
subplot(2,2,1);
for i = 1: 6
    A(i) = plot(prob.time,jnt_init(:,i),'linewidth',2);
    hold on;
end
legend(A,'关节1','关节2','关节3','关节4','关节5','关节6');
grid on;xlim([prob.time(1),prob.time(end)]);
title('(a) 关节位置曲线');
xlabel('时间 T(s)');ylabel('角度 \theta(rad)');

subplot(2,2,2);
for i = 1: 6
    plot(prob.time,vel_init(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([prob.time(1),prob.time(end)]);
title('(b) 关节速度曲线');
xlabel('时间 T(s)');ylabel('速度 v(rad/s)');

subplot(2,2,3);
for i = 1: 6
    plot(prob.time,acc_init(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([prob.time(1),prob.time(end)]);
title('(c) 关节加速度曲线');
xlabel('时间 T(s)');ylabel('加速度 a(rad/s^2)');

subplot(2,2,4);
for i = 1: 6
    plot(prob.time,jerk_init(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([prob.time(1),prob.time(end)]);
title('(d) 关节Jerk曲线');
xlabel('时间 T(s)');ylabel('Jerk j(rad/s^3)');


figure(3)
subplot(1,2,1)
theta=zeros(6,1);
P=zeros(3,length(pos(:,1)));
P_init=zeros(3,length(pos(:,1)));
for i = 1: length(pos(:,1))
    theta(2,1)=-pos(i,2);
    theta(3,1)=-pos(i,3);
    theta(5,1)=-pos(i,5);
    [P(:,i),R]=T12_fkine(theta);
    theta(2,1)=jnt2(i,3);
    theta(3,1)=jnt3(i,3);
    theta(5,1)=jnt5(i,3);
    [P_init(:,i),R]=T12_fkine(theta);
end
plot(prob.time,P_init(1,:),'linewidth',2);
hold on;
plot(prob.time,P(1,:),'linewidth',2);
grid on;xlim([prob.time(1),prob.time(end)]);
legend('原始轨迹','优化后轨迹');
title('(a) 末端位置坐标');
xlabel('时间 T(s)');ylabel('coordinate (m)');
subplot(1,2,2)
plot(prob.time,(P(1,:)-P_init(1,:)).*1000*0.2,'linewidth',2);
grid on;xlim([prob.time(1),prob.time(end)]);
title('(b) 末端误差');
xlabel('时间 T(s)');ylabel('coordinate (mm)');

fid = fopen('pathInfo.txt', 'wt');
for i = 1:length(pos)
    %fprintf(fid,'%f,%f,%f,%f,%f,%f\n',0,-pos(i,2),-pos(i,3),0,-pos(i,5),0);
    fprintf(fid,'%f,%f,%f,%f,%f,%f\n',0,-pos(i,2),-pos(i,3),0,-pos(i,5),0);
end
fclose(fid);
