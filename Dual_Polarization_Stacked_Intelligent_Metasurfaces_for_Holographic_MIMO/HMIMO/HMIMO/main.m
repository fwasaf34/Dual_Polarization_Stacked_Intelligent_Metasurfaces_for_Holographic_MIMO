clc
clear
close all


%%
clearvars -except seed CHmimoIdeal CMimo Loss
rng(1)
Parameter
TX_SIM_Initialize
RX_SIM_Initialize
Channel

[E,Lam,F]=svd(G);
lam=Lam(1:S,1:S)

WaterFilling

%%
for seed=1:100
    rng(seed)
    Phi_Initialize
    TX_SIM
    RX_SIM
    H=GRxSim*G*GTxSim;
    alpha=inv(H(:)'*H(:))*H(:)'*lam(:);
    LossInitialize(seed)=norm(alpha*H-lam,'fro')^2/norm(lam,'fro')^2;
end

%%
LossInitialize;
[~,index]=min(LossInitialize)
rng(index)
Phi_Initialize
TX_SIM
RX_SIM
H=GRxSim*G*GTxSim;
alpha=inv(H(:)'*H(:))*H(:)'*lam(:);
Loss(1)=norm(alpha*H-lam,'fro')^2/norm(lam,'fro')^2
% CHmimo(1)=log2(real(det(eye(S,S)+P/S/noise*H*H')));

%%
count=2;
for round=2:21
    round
    UpdateTxPhi
    UpdateRxPhi
    eata=eata*beta;
    CHmimoCal
end

%%
CHmimoIdeal=sum(log2(1+P/S*diag(lam).^2/noise));
CHmimoIdealWF=sum(log2(1+Ps.*diag(lam).^2/noise));

%plot
% figure
% hold on
% yline(CHmimoIdeal,'-r')
% plot(CHmimo,'-b')
% hold off
% 
figure
hold on
plot(Loss,'-r')
hold off

figure
h=heatmap(abs(alpha*H)); 
% caxis([max(max(abs(H)))*10^-1 max(max(abs(H)))])
h.CellLabelColor = 'none';
% caxis([max(max(abs(alpha*DpH)))*10^-1 max(max(abs(alpha*DpH)))])
caxis([5*10^-8 5*10^-7])
ylabel('Data stream index')
xlabel('Data stream index')

abs(alpha)
CHmimoIdeal(end)
CHmimo(end)
Loss(end)
[CHmimoIdeal CHmimo(end) CHmimoIdealWF CHmimoWF(end)]/P