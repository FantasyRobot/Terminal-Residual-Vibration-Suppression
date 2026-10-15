clear all; clc;
orig = load('原规划振动曲线22.txt');

prob.dt = 0.002;
prob.tf=length(orig(2000:end,5))*prob.dt;
time = 0:prob.dt:prob.tf;
prob.time =time(:,1:end-1);
plot(prob.time,orig(2000:end,5)+9.8);
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
grid on;
xlabel('时间 T(s)');ylabel('加速度 a(g/s^2)');