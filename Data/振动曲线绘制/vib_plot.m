clear all;clc;

vib_init=xlsread("vib_init.xlsx");
vib_opt=xlsread("vib_opt.xlsx");
vib_shape=xlsread("vib_shape.xlsx");
cycle=length(vib_init(:,1));
time=0:0.004:(cycle-1)*0.004;

figure
set(0, 'DefaultAxesFontSize', 20);       % 坐标轴标签和标题字体大小
set(0, 'DefaultTextFontSize', 20);        % 文本对象字体大小
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');  % 使用LaTeX解释器（可选）

plot(time,vib_init(:,5)+9.8,'linewidth',2);
hold on;
plot(time,vib_shape(:,5)+9.8,'linewidth',2);
plot(time,vib_opt(:,5)+9.8,'linewidth',2);
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
hold on;
legend('S型规划振动曲线','输入整形后振动曲线','本文优化后振动曲线');
%legend('Raw vibration','Vibration after input shaping','Vibration after optimization');
grid on;xlim([time(1),time(end)]);
% title('末端Z方向振动加速度曲线');
xlabel('时间 t(s)');ylabel('加速度 a(m/s^2)');
%xlabel('Time t(s)');ylabel('a(m/s^2)');
% zp = BaseZoom();
% zp.plot;