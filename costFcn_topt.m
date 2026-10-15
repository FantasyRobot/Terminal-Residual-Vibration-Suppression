function [ cost ] = costFcn_topt( Xc,Uc,tc,D,nS,nU,ww,scale,prob )
%% cost function of 'pure' time optimal control
cost = scale*2;% only consider trajectory time

%% Xc:npts*nS
Xc = reshape(Xc,nS,numel(Xc)/nS);
Vel = (D*(Xc(1:6,:).')/scale).';
Acc = (D*(Xc(7:12,:).')/scale).';
Jerk = (D*(Acc.')/scale).';

%% close to time optimal, with regulation, improves speed
%% cost0
% rob.r.G = [10     6    42    11     2     2];
% cost = cost ...                
%             + 1*sum(ww(:).'.*( sum(diag(1./prob.rob.r.G.^2)*...
%                      Jerk.^2,1) ))*scale;
%% cost1
% A=[0 6 42 0 6 0];
% cost =sum(ww(:).'.*( sum(diag(1./A.^2)*Jerk.^2,1) ))*scale;
%% cost2
A=[1000 1/50 1/50 1000 1/50 1000];
delta_x= zeros(6,6);
% delta_x = Xc(1:6,:)-prob.x_init(1:6,:);

% cost =cost+...
%                0.3*sum( ww(:).'.*((  sum( diag(A)*(Jerk).^2, 1) + sum( diag(A)*(delta_x).^2,1))))*scale;
 cost =sum(ww(:).'.*( sum(diag(1./A.^2)*Jerk.^2,1) ))*scale;
end
% 
