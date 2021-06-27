set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',15)
%% initialization
%%%% dimension
clear all
% load('ross_rk4_alter_c_matlab2018_downSamp.mat');
load('ross_rk4_alter_c_matlab2018.mat');

Nin=1;
Nout=2;
Nlist = [100 200 400 600 800 1000];
for iN = 1:length(Nlist)
 
    Nres=Nlist(iN)
    alph=0.05;

        
    ru=1;
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

    sigm=1;
    Win=2*sigm*(rand(Nres,Nin)-0.5);

    xi=1;
      clear r

     r(:,1)=rand(Nres,1);

    K1=floor(0.8*K);
    for i=2:K1
        r(:,i)=(1-alph)*r(:,i-1)+alph*tanh(A*r(:,i-1)+Win*u(:,i-1)+xi);

    end

    % optimised minimizers
    clear rbar sbar dS dR s_rc
    bt=0.01;
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


    for i=K1+1:K
        r(:,i)=(1-alph)*r(:,i-1)+alph*tanh(A*r(:,i-1)+Win*u(:,i-1)+xi);
        s_rc(:,i)=Wout*r(:,i)+c;

    end

    RMSE_testy(iN) = sqrt(mean((s_rc(1,K1+1:K)-s(1,K1+1:K)).^2));
    RMSE_testz(iN) = sqrt(mean((s_rc(2,K1+1:K)-s(2,K1+1:K)).^2));


end
%%
plot(Nlist,RMSE_testy,'-x','DisplayName','RMSE_y'); hold on
plot(Nlist,RMSE_testz,'-x','DisplayName','RMSE_z'); hold on
legend
xlabel('Number of nodes');
ylabel('RMSE');
title(['Software RC, \beta = ' num2str(bt) ', \alpha = ' num2str(alph)]);
box on
