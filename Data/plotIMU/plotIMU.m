clear all; clc;
orig = load('original.txt');
cost1 = load('cost1.txt');
jnt = load('oriJnt.txt');

subplot(3,1,1)
plot(orig(:,5)+9.8);
title('原直线规划')
grid on;
subplot(3,1,2)
plot(jnt(:,5)+9.8);
title('原关节规划')
grid on;
subplot(3,1,3)
plot((cost1(:,5)+9.8)*0.08);
title('cost1优化规划')
grid on;