
function[]=magnetic_components_for_Valhalla2()
%XYZ transform to ENR cordinate
clear;
  fpath1='D:\Program Files\MATLAB\R2023b\bin\Crustal field model for Valhalla basin\';
  fpath=[fpath1,'matrix_R\'];
   load([fpath,'Xk'],'Xk');   
K=length(Xk)./3;
d1=Xk(1:K);d2=Xk(K+1:2*K);d3=Xk(2*K+1:3*K); %Mx,My,Mz
di=[d1 d2 d3];
di(:,4)=sqrt(d1.^2+d2.^2+d3.^2); %M_total

A=[];A1=[];A2=[];A3=[];J=0;W=0;dw=2; h=-60;
for j=(270):dw:(350)
    for w=(-20):dw:(40)
       
        J=j;W=w;
        A1=[A1;J,W,h;];
       
    end
end
A=A1; %
di(:,10:11)=A1(:,1:2);

a1=[];
for i=1:length(A1)
    lon=A1(i,1);lat=A1(i,2);
rotation_Ephio_ENR=[-sind(lon),cosd(lon),0;...
    -sind(lat)*cosd(lon),-sind(lat)*sind(lon),cosd(lat);...
    cosd(lat)*cosd(lon),cosd(lat)*sind(lon),sind(lat)]; 
a=rotation_Ephio_ENR*di(i,1:3)'
a1=[a1,a];
end
a1=a1';
a1(:,4)=sqrt(a1(:,1).^2+a1(:,2).^2+a1(:,3).^2);% Totle field

di(:,5:7)=a1(:,1:3);%%ENR components
di(:,8)=asind(di(:,7)./di(:,4)); %Inclination
di(:,10)=sqrt(di(:,5).^2+di(:,6).^2);%Horizontal components
for i=1:length(di)
    if di(i,1)>=0
di(i,9)=acosd(di(i,6)./di(i,10)); % Declination
    else
di(i,9)=360-acosd(di(i,6)./di(i,10));      
    end
end
di1=di;
save([fpath1,'1271_dipoles\di1.mat'],'di1');


figure(4)
subplot('Position',[0.07,0.575,0.43,0.37])
h0=scatter(A1(:,1),A1(:,2),150,di(:,4),'filled');grid on;
c1=colorbar;
c1.Title.String='Dipole Moment (Am^{2})';
c1.Title.FontSize=8;
xlim([267,352]);
ylim([-22,42]);
xlabel('Longitude (бу)','FontSize',8);   
ylabel('Latitude (бу)','FontSize',8); 
set(gca,'FontSize',8);



subplot('Position',[0.57,0.575,0.43,0.37])
h3=scatter(A1(:,1),A1(:,2),150,di(:,9),'filled');grid on;
c3=colorbar;
c3.Title.String='Declination (бу)';
c3.Title.FontSize=8;
xlim([267,352]);
ylim([-22,42]);
xlabel('Longitude (бу)','FontSize',8);   
ylabel('Latitude (бу)','FontSize',8); 
set(gca,'FontSize',8);

subplot('Position',[0.07,0.095,0.43,0.37])
h1=scatter(A1(:,1),A1(:,2),150,di(:,8),'filled');grid on;
c2=colorbar;
c2.Title.String='Inclination (бу)';
c2.Title.FontSize=8;
xlim([267,352]);
ylim([-22,42]);
xlabel('Longitude (бу)','FontSize',8);   
ylabel('Latitude (бу)','FontSize',8); 
set(gca,'FontSize',8);

