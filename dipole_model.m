function[x0,y0,z0]=dipole_model(Jmin,Jmax,Wmin,Wmax,dw,h)



A=[];A1=[];A2=[];A3=[];J=0;W=0;
for j=(Jmin):dw:(Jmax)
    for w=(Wmin):dw:(Wmax)
       
        J=j;W=w;
        A1=[A1;J,W,h;];
       
    end
end
A=A1;
%[£¨J£¬W£¬H£©->(Gx,Gy,Gz)]
J=A(:,1);
W=A(:,2);
H=A(:,3);
n=length(H);
D=ones(n,1);
D=2409*D;
        r=D+H;
        x0=1000*r.*cos(2*pi*W/360).*cos(2*pi*J/360);
        y0=1000*r.*cos(2*pi*W/360).*sin(2*pi*J/360);
        z0=1000*r.*sin(2*pi*W/360);

        