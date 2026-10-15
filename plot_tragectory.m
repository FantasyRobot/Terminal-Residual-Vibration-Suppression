pos = load("pathInfo.txt");
pos_time = zeros(length(pos), 1);
vel_jnt = zeros(length(pos), 6);
acc_jnt = zeros(length(pos), 6);
jerk_jnt = zeros(length(pos), 6);


jnt2 = load("jnt2.txt");
jnt3 = load("jnt3.txt");
jnt5 = load("jnt5.txt");
data_length = length(pos_time);
jnt_init(:,1) = zeros(data_length,1);
jnt_init(:,2) = -jnt2(1:data_length,3);
jnt_init(:,3) = -jnt3(1:data_length,3);
jnt_init(:,4) = zeros(data_length,1);
jnt_init(:,5) = -jnt5(1:data_length,3);
jnt_init(:,6) = zeros(data_length,1);
vel_init=zeros(data_length,6);
acc_init=zeros(data_length,6);
jerk_init=zeros(data_length,6);

for i = 1: length(pos_time)
    pos_time(i,1)=i*0.004;
    hold on;
end

figure (1)
subplot(2,2,1);
for i = 1: 6
    pos(:,i) = -pos(:,i);
    A(i) = plot(pos_time, pos(:,i),'linewidth',2);
    hold on;
end
legend(A,'Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6');
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('\theta(бу)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(a) Position');

subplot(2,2,2);
for i = 1: length(vel_jnt)-1
    for j = 1: 6
        vel_jnt(i,j) = pos(i+1,j) - pos(i,j);
    end
end
for i = 1: 6
    plot(pos_time, vel_jnt(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('v(бу/s)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(b) Velocity');

subplot(2,2,3);
for i = 1: length(vel_jnt)-1
    for j = 1: 6
        acc_jnt(i,j) = vel_jnt(i+1,j) - vel_jnt(i,j);
    end
end
for i = 1: 6
    plot(pos_time, acc_jnt(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('a(бу/s^2)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(c) Acceleration');

subplot(2,2,4);
for i = 1: length(acc_jnt)-1
    for j = 1: 6
        jerk_jnt(i,j) = acc_jnt(i+1,j) - acc_jnt(i,j);
    end
end
for i = 1: 6
    plot(pos_time, jerk_jnt(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('a(бу/s^3)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(c) Jerk');

figure (2)
subplot(2,2,1);
for i = 1: 6
    A(i) = plot(pos_time, jnt_init(:,i),'linewidth',2);
    hold on;
end
legend(A,'Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6');
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('\theta(бу)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(a) Position');

subplot(2,2,2);
for i = 1: length(jnt_init)-1
    for j = 1: 6
        vel_init(i,j) = jnt_init(i+1,j) - jnt_init(i,j);
    end
end
for i = 1: 6
    plot(pos_time, vel_init(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('v(бу/s)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(b) Velocity');

subplot(2,2,3);
for i = 1: length(vel_init)-1
    for j = 1: 6
        acc_init(i,j) = vel_init(i+1,j) - vel_init(i,j);
    end
end
for i = 1: 6
    plot(pos_time, acc_init(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('a(бу/s^2)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(c) Acceleration');

subplot(2,2,4);
for i = 1: length(acc_init)-1
    for j = 1: 6
        jerk_init(i,j) = acc_init(i+1,j) - acc_init(i,j);
    end
end
for i = 1: 6
    plot(pos_time, jerk_init(:,i),'linewidth',2);
    hold on;
end
grid on;xlim([0,3.8]);
xlabel('Time T(s)');ylabel('a(бу/s^3)');
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
title('(c) Jerk');