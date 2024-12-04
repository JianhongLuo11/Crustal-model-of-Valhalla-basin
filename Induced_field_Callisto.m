%Magnetic field induced from the subsurface ocean of Callisto

% rx,ry,rz:spacecraft positions in CphiO cordinate 
% Hdipole:Equatorial Field Magnitude in nT
% theta_z,theta_y: Angle From spin axis of Callisto and Angle From Jupiter Facing Meridian
% Longitude:Spacecraft Longitude

% The induced fields for the C3 flyby are shown as follows£º
% induced_field_C3=Induced_field_Callisto(data_C3_induce(:,1),data_C3_induce(:,2),data_C3_induce(:,3),15,90,7,data_C3_induce(:,4))  

function [Bx,By,Bz,Bt]=Induced_field_Callisto(rx,ry,rz,Hdipole,theta_z,theta_y,Longitude)
[x1,y1,z1,t1]=Induced_field_Z(rx,ry,rz,Hdipole,theta_z,theta_y);
[x2,y2,z2,t2]=Induced_field_Y(rx,ry,rz,Hdipole,theta_z,theta_y);
[x3,y3,z3,t3]=Induced_field_X(rx,ry,rz,Hdipole,theta_z,theta_y);
z_conponent=[x1,y1,z1];% magnetic field components magnetised along -Z axis
y_conponent=[x2,y2,z2];% -Y axis
x_conponent=[x3,y3,z3];% -X axis
if theta_z>=90 %-z
    z_conponent=z_conponent;
elseif theta_z<90
    z_conponent=-z_conponent;%z
end
if theta_y>=0 && theta_y<90 %x,y
    y_conponent=-y_conponent;
    x_conponent=-x_conponent;
elseif theta_y>=90 && theta_y<180 %x,-y
    y_conponent=y_conponent;
    x_conponent=-x_conponent;
elseif theta_y>=180 && theta_y<270 %-x,-y
    y_conponent=y_conponent;
    x_conponent=x_conponent;
elseif theta_y>=270 && theta_y<360 %-x,y
    y_conponent=-y_conponent;
    x_conponent=x_conponent;
end
Bx=z_conponent(:,1)+y_conponent(:,1)+x_conponent(:,1);% X components arised from the sum of three dipoles
By=z_conponent(:,2)+y_conponent(:,2)+x_conponent(:,2);
Bz=z_conponent(:,3)+y_conponent(:,3)+x_conponent(:,3);
Bt=sqrt(Bx.^2+By.^2+Bz.^2);
plot(Longitude(:),Bt(:),'ko');hold on;grid on;
xlim([180,360]);
xlabel('Longitude(¡ã)');
ylabel('B(nT)');
plot(Longitude(:),Bx(:),'bo');hold on;
plot(Longitude(:),By(:),'mo');hold on
plot(Longitude(:),Bz(:),'go');
legend('B_t_o_t_a_l','B_x','B_y','B_z')
end


%subfunctions
% dipole magnetized along -Z axis
function [Hx,Hy,Hz,Ht] = Induced_field_Z(rx,ry,rz,Hdipole,theta_z,theta_y)
% %-Z
Hx=[];
Hy=[];
Hz=[];
% Ht=[];
i=1;
while i<=length(rx)
    Rx=rx(i);
    Ry=ry(i);
    Rz=rz(i);
AB=sqrt((Rx).^2+(Ry).^2+(Rz).^2);
Lat=asind(Rz./AB);%Latitude£¬-Z
I=atand(2*tand(Lat));%
H=abs(Hdipole*cosd(theta_z)*sqrt(3*(sind(Lat)).^2+1)*((1./AB)^3));%
Hh=H*cosd(I);%Horizental components
Hv=H*sind(I);%Vertical components
Hvz=-Hv*sind(Lat);%
Hvh=Hv*cosd(Lat);%
Hhz=abs(Hh*cosd(Lat));% 
Hhh=Hh*sind(Lat);%

theta=atand(Ry./Rx);
if Rz>=0 && Rx>=0 && Ry>=0 %1
    Hvx=-abs(Hvh*cosd(theta));
    Hvy=-abs(Hvh*sind(theta));
    Hhx=-abs(Hhh*cosd(theta));
    Hhy=-abs(Hhh*sind(theta));
elseif Rz>=0 && Rx<0 && Ry>=0 %2
    Hvx=abs(Hvh*cosd(theta));
    Hvy=-abs(Hvh*sind(theta));
    Hhx=abs(Hhh*cosd(theta));
    Hhy=-abs(Hhh*sind(theta));
elseif Rz>=0 && Rx<0 && Ry<0 %3
    Hvx=abs(Hvh*cosd(theta));
    Hvy=abs(Hvh*sind(theta));
    Hhx=abs(Hhh*cosd(theta));
    Hhy=abs(Hhh*sind(theta));
elseif Rz>=0 && Rx>=0 && Ry<0 %4
    Hvx=-abs(Hvh*cosd(theta));
    Hvy=abs(Hvh*sind(theta));
    Hhx=-abs(Hhh*cosd(theta));
    Hhy=abs(Hhh*sind(theta));
elseif Rz<0 && Rx>=0 && Ry>=0 %5
    Hvx=abs(Hvh*cosd(theta));
    Hvy=abs(Hvh*sind(theta));
    Hhx=abs(Hhh*cosd(theta));
    Hhy=abs(Hhh*sind(theta));
elseif Rz<0 && Rx<0 && Ry>=0 %6
    Hvx=-abs(Hvh*cosd(theta));
    Hvy=abs(Hvh*sind(theta));
    Hhx=-abs(Hhh*cosd(theta));
    Hhy=abs(Hhh*sind(theta));
elseif Rz<0 && Rx<0 && Ry<0 %7
    Hvx=-abs(Hvh*cosd(theta));
    Hvy=-abs(Hvh*sind(theta));
    Hhx=-abs(Hhh*cosd(theta));
    Hhy=-abs(Hhh*sind(theta));
elseif Rz<0 && Rx>=0 && Ry<0 %8
    Hvx=abs(Hvh*cosd(theta));
    Hvy=-abs(Hvh*sind(theta));
    Hhx=abs(Hhh*cosd(theta));
    Hhy=-abs(Hhh*sind(theta));
    
end
Hx=[Hx;Hvx+Hhx];
Hy=[Hy;Hvy+Hhy];
Hz=[Hz;Hvz+Hhz];
i=i+1;
end    
Ht=sqrt(Hx.^2+Hy.^2+Hz.^2);

% plot(Longitude(:),Ht(:),'ko-');hold on;grid on;
% xlim([200,360]);
% xlabel('Longitude(¡ã)');
% ylabel('magnetic field(nT)');
% plot(Longitude(:),Hx(:),'bo-');hold on;
% plot(Longitude(:),Hy(:),'mo-');hold on
% plot(Longitude(:),Hz(:),'go-');
% legend('B_t_o_t_a_l','B_x','B_y','B_z')
end

% -Y axis
function [Hx,Hy,Hz,Ht] = Induced_field_Y(rx,ry,rz,Hdipole,theta_z,theta_y)
%-Y
Hx=[];
Hy=[];
Hz=[];
% Ht=[];
i=1;
while i<=length(rx)
    Rx=rx(i);
    Ry=ry(i);
    Rz=rz(i);
AB=sqrt((Rx).^2+(Ry).^2+(Rz).^2);
Lat=asind(Ry./AB);
I=atand(2*tand(Lat));
H=abs(Hdipole*sind(theta_z)*cosd(theta_y)*sqrt(3*(sind(Lat))^2+1)*((1./AB)^3));%
Hh=H*cosd(I);
Hv=H*sind(I);

Hvy=-abs(Hv*sind(Lat));
Hvv=Hv*cosd(Lat);
Hhy=abs(Hh*cosd(Lat));
Hhh=Hh*sind(Lat);
theta=atand(Rz/Rx);

if Ry>=0 && Rz>=0 && Rx>=0
Hvvz=-abs(Hvv*sind(theta));
Hvvx=-abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhx=-abs(Hhh*cosd(theta));
elseif Ry>=0 && Rz>=0 && Rx<0 
Hvvz=-abs(Hvv*sind(theta));
Hvvx=abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhx=abs(Hhh*cosd(theta));
elseif Ry>=0 && Rz<0 && Rx<0 
Hvvz=abs(Hvv*sind(theta));
Hvvx=abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhx=abs(Hhh*cosd(theta));
elseif Ry>=0 && Rz<0 && Rx>=0 
Hvvz=abs(Hvv*sind(theta));
Hvvx=-abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhx=-abs(Hhh*cosd(theta));
  elseif  Ry<0 && Rz>=0 && Rx>=0 
Hvvz=abs(Hvv*sind(theta));
Hvvx=abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhx=abs(Hhh*cosd(theta));
  elseif Ry<0 && Rz>=0 && Rx<0 
Hvvz=abs(Hvv*sind(theta));
Hvvx=-abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhx=-abs(Hhh*cosd(theta));
  elseif  Ry<0 && Rz<0 && Rx<0 
Hvvz=-abs(Hvv*sind(theta));
Hvvx=-abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhx=-abs(Hhh*cosd(theta));
  elseif Ry<0 && Rz<0 && Rx>=0 
Hvvz=-abs(Hvv*sind(theta));
Hvvx=abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhx=abs(Hhh*cosd(theta));
end
Hx=[Hx;Hvvx+Hhhx];
Hy=[Hy;Hvy+Hhy];
Hz=[Hz;Hvvz+Hhhz];
i=i+1;
end
Ht=sqrt(Hx.^2+Hy.^2+Hz.^2);
% plot(Longitude,Ht,'ko-');hold on;grid on;
% xlim([200,360]);
% xlabel('Longitude(¡ã)');
% ylabel('magnetic field(nT)');
% plot(Longitude,Hx,'bo-');hold on;
% plot(Longitude,Hy,'mo-');hold on
% plot(Longitude,Hz,'go-');
% legend('B_t_o_t_a_l','B_x','B_y','B_z')
end

% -X axis
function [Hx,Hy,Hz,Ht] = Induced_field_X(rx,ry,rz,Hdipole,theta_z,theta_y)
Hx=[];
Hy=[];
Hz=[];
% Ht=[];
i=1;
while i<=length(rx)
    Rx=rx(i);
    Ry=ry(i);
    Rz=rz(i);
AB=sqrt((Rx).^2+(Ry).^2+(Rz).^2);
Lat=asind(Rx./AB);%
I=atand(2*tand(Lat));
H=abs(Hdipole*sind(theta_z)*sind(theta_y)*sqrt(3*(sind(Lat))^2+1)*((1./AB)^3));%
Hh=H*cosd(I);%
Hv=H*sind(I);%

Hvx=-abs(Hv*sind(Lat));
Hvv=Hv*cosd(Lat);
Hhx=abs(Hh*cosd(Lat));
Hhh=Hh*sind(Lat);
theta=atand(Rz/Ry);

if Rx>=0 && Ry>=0 && Rz>=0  %+x
Hvvz=-abs(Hvv*sind(theta));
Hvvy=-abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhy=-abs(Hhh*cosd(theta));   
elseif Rx>=0 && Ry<0 && Rz>=0 
Hvvz=-abs(Hvv*sind(theta));
Hvvy=abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhy=abs(Hhh*cosd(theta));    
elseif Rx>=0 && Ry<0 && Rz<0  
Hvvz=abs(Hvv*sind(theta));
Hvvy=abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhy=abs(Hhh*cosd(theta));    
elseif Rx>=0 && Ry>=0 && Rz<0  
Hvvz=abs(Hvv*sind(theta));
Hvvy=-abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhy=-abs(Hhh*cosd(theta)); 
elseif Rx<0 && Ry>=0 && Rz>=0  
Hvvz=abs(Hvv*sind(theta));
Hvvy=abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhy=abs(Hhh*cosd(theta));    
elseif Rx<0 && Ry<0 && Rz>=0  
Hvvz=abs(Hvv*sind(theta));
Hvvy=-abs(Hvv*cosd(theta));
Hhhz=abs(Hhh*sind(theta));
Hhhy=-abs(Hhh*cosd(theta));    
elseif Rx<0 && Ry<0 && Rz<0  
Hvvz=-abs(Hvv*sind(theta));
Hvvy=-abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhy=-abs(Hhh*cosd(theta));    
elseif Rx<0 && Ry>=0 && Rz<0   
Hvvz=-abs(Hvv*sind(theta));
Hvvy=abs(Hvv*cosd(theta));
Hhhz=-abs(Hhh*sind(theta));
Hhhy=abs(Hhh*cosd(theta));     
end
Hx=[Hx;Hvx+Hhx];
Hy=[Hy;Hvvy+Hhhy];
Hz=[Hz;Hvvz+Hhhz];
i=i+1;
end
Ht=sqrt(Hx.^2+Hy.^2+Hz.^2);
% plot(Longitude,Ht,'ko-');hold on;grid on;
% xlim([200,360]);
% xlabel('Longitude(¡ã)');
% ylabel('magnetic field(nT)');
% plot(Longitude,Hx,'bo-');hold on;
% plot(Longitude,Hy,'mo-');hold on
% plot(Longitude,Hz,'go-');
% legend('B_t_o_t_a_l','B_x','B_y','B_z')
end



