function [ lb,ub,bnds ] = bounds( npts,init,target,prob )
% generate control lower bounds and upper bounds

rob = prob.rob;

% bounds for velocity, torque, and torque rate
posinf=rob.r.Qlim(1,:).'; 
possup=rob.r.Qlim(2,:).';

velbnd=rob.r.QDlim(2,:).';

% trqbnd=0.8*[1396.5,1402.3,382.7,45.2,44.6,32.5].';
% dtrqbnd=30*trqbnd;
trqbnd=0.1*[1396.5,1402.3,382.7,45.2,44.6,32.5].';
dtrqbnd=0.1*trqbnd;

bnds = [posinf,possup,velbnd,trqbnd,dtrqbnd];

limit=1;% of limit of every bound
interpLimit=0.7; % for test case
% interpLimit=0.68;% considering interpolation using 20 pnts, addition limit
% interpLimit=0.7;% considering interpolation using 12 pnts, addition limit
conlb = -[-posinf;velbnd;trqbnd;dtrqbnd]*limit*interpLimit;
conub =  [ possup;velbnd;trqbnd;dtrqbnd]*limit*interpLimit;

tlimit = 100;
lb = [kron(ones(npts,1),conlb);...
    init(:);
    zeros(6,1);...
    zeros(6,1);...
    zeros(6,1);...
    target(:);...
    zeros(6,1);...
    zeros(6,1);...
    zeros(6,1);...
    0];

ub = [kron(ones(npts,1),conub);...
    init(:);
    zeros(6,1);...
    zeros(6,1);...
    zeros(6,1);...
    target(:);...
    zeros(6,1);...
    zeros(6,1);...
    zeros(6,1);...
    tlimit];

% self collision
lb=[lb;0.0*ones(npts*size(prob.selfmap,2),1)];
ub=[ub;10000*ones(npts*size(prob.selfmap,2),1)];
% wall collision
wlb=[];
wub=[];
for k=1:size(prob.wallmap,2)
    Lind=prob.wallmap(k).ind(1);
    Bind=prob.wallmap(k).ind(2);
    r=prob.r{Lind}(Bind);
    wlb=[wlb;r+0.05];
    wub=[wub;10000];
end
lb=[lb;kron(ones(npts,1),wlb)];
ub=[ub;kron(ones(npts,1),wub)];
% obstacle collision
lb=[lb;0.003*ones(npts*size(prob.obsmap,2),1)];
ub=[ub;10000*ones(npts*size(prob.obsmap,2),1)];

end

