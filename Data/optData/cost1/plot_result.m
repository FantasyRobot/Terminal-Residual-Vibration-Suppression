jntTorqure=load('torqueValue.txt');
plot(jntTorqure);
legend('关节1','关节2','关节3','关节4','关节5','关节6');
jntRa=load('jointValue.txt');
plot(jntRa);
legend('关节1','关节2','关节3','关节4','关节5','关节6');
for j =1:6
    for i=1:length(jntRa)-1
        vel(i,j)=jntRa(i+1,j)-jntRa(i,j);
    end
end
plot(vel);
legend('关节1','关节2','关节3','关节4','关节5','关节6');