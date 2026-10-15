clear;
clc;
% f = chebfun('abs(x-.3)');
% % 函数长度：
% length(f)

% f = chebfun('abs(x-.3)','splitting', 'on')
% 三个函数 x^{2},1,4-x ，以及一个端点向量，表示第一个函数应用于 [-1,1] ，第二个函数应用于 [1,2] ，第三个函数应用于 [2,4] ：
f = chebfun({@(x) x^2, 1, @(x) 4-x},[-1 1 2 4]);
plot(f)
