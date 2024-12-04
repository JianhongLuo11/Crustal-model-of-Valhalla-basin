function[]=main_gongetidu_3D_2()
% clear;clc;close all
%calculating the Xk(Mx,My,Mz)
fpath1='D:\Program Files\MATLAB\R2023b\bin\Crustal field model for Valhalla basin\';

G=[];
for kk=1:1:3
    %
    switch kk
        case 1
            fpath=[fpath1,'matrix_R\x\'];
        case 2
            fpath=[fpath1,'matrix_R\y\'];
        case 3
            fpath=[fpath1,'matrix_R\z\'];
    end
    load([fpath1,'matrix_R\b'],'b');
 
    load([fpath,'Gx'],'Gx');
    load([fpath,'Gy'],'Gy');
    load([fpath,'Gz'],'Gz');
    G1=[Gx;Gy;Gz]; %%%%
    A1=size(G1); %%
    G=[G,G1];
end
    1
    A2=size(G);
    
    Q=G'*G;
    B=G'*b;

    2
    
    
    a=length(G(1,:));
    n=length(b(:,1));
%     n_c=length(b_c(:,1))
    m0=1;
    X0=ones(a,1)*m0; %%initial value of m

    k=[];
    aa=100; %% 100%
    x0=X0; 
    r0=B-Q*x0; 
    d0=r0;
    k=0;
    %rk is the error of the kth iteration, and dk is the conjugate vector we require; k is the number of iterations
    rk=r0;dk=d0;xk=x0;s=[];xk1=x0;
    3


    while aa>0.5 %%0.5
        
        a1=(rk'*rk)/(dk'*Q*dk); %A5
        xk1=xk+a1*dk;
        rk1=B-Q*xk1;
        bx=G*xk1;
        www=b-bx;
        
       s(k+1)=sqrt((b-bx)'*(b-bx)/(n)); %%

  
%         bx_c=G_c*xk1;
%         s1(k+1)=sqrt((b_c-bx_c)'*(b_c-bx_c)/n_c);
        
        if k>0
            aa=(abs(s(k))-abs(s(k+1)))/abs(s(k))*100
        end
        b1=(rk1'*rk1)/(rk'*rk);
        dk1=rk1+b1*dk;
        rk=rk1;
        dk=dk1;
        k=k+1;
        xk=xk1;
        
    end
    4
    Xk=xk;
    save([fpath1,'matrix_R\Xk.mat'],'Xk');


    save([fpath1,'matrix_R\s.mat'],'s');

    smin=find(s==min(s));
%     ee=min(s);
% eee=length(s);
    x=1:1:k;
    figure(kk)
    scatter(x,s,300,'.');grid on;
% 
% 
    xlabel('Numbers of Iteration','fontsize',12);
    ylabel('RMS (nT)','fontsize',12);
%     clear;
