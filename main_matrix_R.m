function[]=main_matrix_R()
%%%b=Gm+v this function calculated the G value

%2023.12.15
% clear;clc;

load data_R_fun.mat;  % fit data
% column 1-9 represent longitude, latitude, height, magnetic three components in PC cordinate (planetocentric coordinate in Li et al., 2020),
% positions in x,y,z in PC cordinate.
% longitude, latitude, height, GBx,GBy,GBz,Gx,Gy,Gz,respectively (see Li et al., 2020).

%longitude and latitude range of the dipoles
        Jmin=270;Jmax=350;Wmin=-20;Wmax=40;
        dw=2; h=-60;%steps and burid depth

%fpath='D:\model-Lxz-2020';
%milename='';
%load([fpath,milename]);
% fpath1='D:\matlab2018a\bin\callisto\';
fpath1='D:\Program Files\MATLAB\R2023b\bin\Crustal field model for Valhalla basin\';



%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+
%Gx=b£»
%Next, construct the b-matrix
Bx=data_R_fun(:,4);
By=data_R_fun(:,5);
Bz=data_R_fun(:,6);
b=[Bx;By;Bz];
save([fpath1,'matrix_R\b.mat'],'b');
b=[];


%calculating the x by  the conjugate gradient technique
[x0 y0 z0]=dipole_model(Jmin,Jmax,Wmin,Wmax,dw,h); 
D=[x0 y0 z0];  %%dipole locations
a1=length(D(:,1))
save([fpath1,'matrix_R\dipole_location.mat'],'D');
D=[];


%construct the G-matrix
x=data_R_fun(:,7);
y=data_R_fun(:,8);
z=data_R_fun(:,9);
R=sqrt(x.^2+y.^2+z.^2);


m=length(x)
n=length(x0)


x0=x0';y0=y0';z0=z0';
XX=bsxfun(@minus,x,x0);
YY=bsxfun(@minus,y,y0);
ZZ=bsxfun(@minus,z,z0);

 r=sqrt(XX.^2+YY.^2+ZZ.^2);
 r3=1./(r.^3);r5=1./(r.^5);

for kk=1:1:3
    %
    switch kk
        case 1
            Gx=-100.*r3+300.*XX.*XX.*r5;
          
            Gy=300.*XX.*YY.*r5;
        
            Gz=300.*XX.*ZZ.*r5;
          
            fpath2=[fpath1,'matrix_R\x\']; %%
        case 2
            Gx=300.*XX.*YY.*r5;
        
            Gy=-100.*r3+300.*YY.*YY.*r5;
         
            Gz=300.*YY.*ZZ.*r5;
        
            fpath2=[fpath1,'matrix_R\y\']; %%
        case 3
            Gx=300.*XX.*ZZ.*r5;
          
            Gy=300.*YY.*ZZ.*r5;
         
            Gz=-100.*r3+300.*ZZ.*ZZ.*r5;
         
            fpath2=[fpath1,'matrix_R\z\']; %%
    end
    save([fpath2,'Gx'],'Gx');
    save([fpath2,'Gy'],'Gy');
    save([fpath2,'Gz'],'Gz');
    Gx=[];Gy=[];Gz=[];
end    


