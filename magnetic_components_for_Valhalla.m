function[]=magnetic_components_for_Valhalla()
%this function calculates the magnetic conponents in PC corninate
% X in PC cordinate equal to Y in CphiO cordinate; Y in PC equal to -X in CphiO;
%Z in PC equal to Z in CphiO.
clear;clc;close all;
        
        Jmin=270;Jmax=350;Wmin=-20;Wmax=40;  
        dw=2; 
        h=-60;
        fpath1='D:\Program Files\MATLAB\R2023b\bin\Crustal field model for Valhalla basin\';
         
    fpath=[fpath1,'matrix_R\'];

    load([fpath,'Xk'],'Xk');
    load([fpath1,'matrix_R\dipole_location'],'D');
    load data_R_fun.mat
    load data_R_fun2.mat
G1=[];G2=[]; 
C=[];
C=data_R_fun2(:,7:9);  

for kk=1:1:3 


    x=C(:,1);
    y=C(:,2);
    z=C(:,3);
    R=sqrt(x.^2+y.^2+z.^2);
    x0=D(:,1); 
    y0=D(:,2);
    z0=D(:,3);
    
    
    m=length(x) %%N
    n=length(x0) %%K
    G1=zeros(m,n);
    
    x0=x0';y0=y0';z0=z0';
    XX=bsxfun(@minus,x,x0);
    YY=bsxfun(@minus,y,y0);
    ZZ=bsxfun(@minus,z,z0);
    r=sqrt(XX.^2+YY.^2+ZZ.^2);
    r3=1./(r.^3);r5=1./(r.^5);   
    switch kk %%
        case 1
            Gx=-100.*r3+300.*XX.*XX.*r5;
            3.11
            Gy=300.*XX.*YY.*r5;
            3.12
            Gz=300.*XX.*ZZ.*r5;
            3.13
        case 2
            Gx=300.*XX.*YY.*r5;
            3.21
            Gy=-100.*r3+300.*YY.*YY.*r5;
            3.22
            Gz=300.*YY.*ZZ.*r5;
            3.23
        case 3
           Gx=300.*XX.*ZZ.*r5;
            3.31
            Gy=300.*YY.*ZZ.*r5;
            3.32
            Gz=-100.*r3+300.*ZZ.*ZZ.*r5;
            3.33
    end
    G1=[Gx Gy Gz];
    G2=[G2;G1];
   
end
   bx=G2*Xk;
   bxyz=[]; 
bxyz(1:length(C),1)=bx(1:length(C));
bxyz(1:length(C),2)=bx(length(C)+1:2*length(C));
bxyz(1:length(C),3)=bx(2*length(C)+1:3*length(C));
bxyz(1:length(C),4)=sqrt(bxyz(:,1).^2+bxyz(:,2).^2+bxyz(:,3).^2);
bxyz1=bxyz;
save([fpath1,'1271_dipoles\bxyz1.mat'],'bxyz1');


