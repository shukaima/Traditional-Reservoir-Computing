set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',15)
%% initialization
%%%% dimension
Nin=1;
Nout=2;
Nres=100;

ru=0.01;
sparserat=0.2;
D=sparserat*Nres;
A = 2*(rand(Nres,Nres)-0.5);
for i=1:Nres
    for j=1:Nres
        x=rand;
        if x>sparserat
            A(i,j)=0;
        end
    end
end

[e2]=eig(A);
b=abs(e2);
A=A.*(ru/max(b));
clear e2 b

sigm=0.01;
Win=2*sigm*(rand(Nres,Nin)-0.5);

xi=0;

alph=0.5;


  clear r

 r(:,1)=rand(Nres,1);

K1=0.8*K;
for i=2:K1
    r(:,i)=(1-alph)*r(:,i-1)+alph*relu_og(A*r(:,i-1)+Win*u(:,i-1)+xi);

end

% optimised minimizers
clear rbar sbar dS dR s_rc
bt=0.1;
btt=bt*eye(Nres);
rbar=mean(r,2);
sbar=mean(s,2);
dS=s-sbar; dR=r-rbar;
Wout=dS(:,1:K1)*dR.'/(dR*dR.'+btt);
c=-(Wout*rbar-sbar);



% for the previous steps
for i=1:K1
    s_rc(:,i)=Wout*r(:,i)+c;

end
% %%
% mk=2;
% RMSE_train = sqrt(mean((s_rc(mk,1:K1)-s(mk,1:K1)).^2))
% 
% h=figure
% hold on
% fig2=plot(1:K1,s(mk,1:K1),'Displayname','from data');
% fig1=plot(1:K1,s_rc(mk,1:K1),':','Displayname','from RC');
% legend([fig1,fig2]);
% %title(['Software RC, train, N =' num2str(Nres) ', \beta = ' num2str(bt) ', RMSE = ' num2str(RMSE_train)]);
% title(['Software RC, train, N =' num2str(Nres) ', \alpha = ' num2str(alph) ', RMSE = ' num2str(RMSE_train)]);
% xlabel('time');
% box on
% ylabel('normalized output and data');
% hold off
% 
% %% let's do the prediction

for i=K1+1:K
    r(:,i)=(1-alph)*r(:,i-1)+alph*tanh(A*r(:,i-1)+Win*u(:,i-1)+xi);
    s_rc(:,i)=Wout*r(:,i)+c;

end
%
mk=1;
RMSE_test = sqrt(mean((s_rc(mk,K1+1:K)-s(mk,K1+1:K)).^2))

h=figure
hold on
fig2=plot(K1+1:K,s(mk,K1+1:K),'Displayname','from data');
fig1=plot(K1+1:K,s_rc(mk,K1+1:K),':','Displayname','from RC');
legend([fig1,fig2]);
title(['Software RC, test, N =' num2str(Nres) ', \beta = ' num2str(bt) ', RMSE = ' num2str(RMSE_test)]);
xlabel('time');
box on
ylabel('normalized output and data');
hold off

 
 