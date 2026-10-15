function [ defects ] = conEq( Xc,Uc,D,nS,nU,scale,P,prob )

% equality constraint
Xc = reshape(Xc,nS,numel(Xc)/nS);
Uc = reshape(Uc,nU,numel(Uc)/nU);

Vel = (D*(Xc(1:6,:).')/scale).';
Acc = (D*(Xc(7:12,:).')/scale).';
Jerk = (D*(Acc.')/scale).'; 
% D*X=(t_f/2)*F
defects = (D*(Xc.')/scale).';
defects = defects(:);

defects = defects - robDyn(Xc,Uc,prob.rob);
end

