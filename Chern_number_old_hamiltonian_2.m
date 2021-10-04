close all
clear 
clc
%%%%%%%%%%%%%%%
L=1;
thetaa=0.3*pi;
thetab=0.01*pi;

Bg1=-0.5*pi/L;
Bg3=-Bg1;
Bg2=1*pi/L;

Bg_test=Bg2;

ka=thetaa/(L/4);
kb=thetab/(L/4);


stepkx=0.05*pi;
stepky=0.05*pi;

dkx=0.001*stepkx;
dky=0.001*stepky;


k_start=-0.5;
k_end=0.5;

kx=(k_start*pi:stepkx:k_end*pi);
ky=(k_start*pi:stepky:k_end*pi);




At=zeros(length(kx),length(ky));

jjx=1;
jjy=1;

   
for ky=(k_start*pi:stepkx:k_end*pi)
    jjx=1;
    for kx=(k_start*pi:stepkx:k_end*pi)
        
H_1=-[0, ka.*exp(1i*kx),0, 0
    ka.*exp(-1i*kx), 0, 0, 0
     0, 0, 0, kb.*exp(1i*kx)
     0, 0, kb.*exp(-1i*kx), 0];
H_2=-[0, 0,ka.*exp(1i*ky), 0
    0, 0, 0, kb.*exp(1i*ky)
     ka.*exp(-1i*ky), 0, 0, 0 
     0, kb.*exp(-1i*ky), 0, 0];
H_3=-[0, ka.*exp(-1i*kx),0, 0
    ka.*exp(1i*kx), 0, 0, 0
     0, 0, 0, kb.*exp(-1i*kx)
     0, 0, kb.*exp(1i*kx), 0];
H_4=-[0, 0,ka.*exp(-1i*ky), 0
    0, 0, 0, kb.*exp(-1i*ky)
     ka.*exp(1i*ky), 0, 0, 0 
     0, kb.*exp(1i*ky), 0, 0];
 
 U1=expm(-1i.*H_1.*L./4);
 U2=expm(-1i.*H_2.*L./4);
 U3=expm(-1i.*H_3.*L./4);
 U4=expm(-1i.*H_4.*L./4);
 U_F=U4*U3*U2*U1;  
 
 
 

   
   
   [V,D]=eig(U_F);
   Q=V;

   En(1)=-log(D(1,1))./1i;
   En(2)=-log(D(2,2))./1i;
   En(3)=-log(D(3,3))./1i;
   En(4)=-log(D(4,4))./1i;
   
   for jj=1:3
   for ii=jj+1:4
       if real(En(jj))>real(En(ii))
           holder_e=En(jj);
           En(jj)=En(ii);
           En(ii)=holder_e;
           holder_v=Q(:,jj);
           Q(:,jj)=Q(:,ii);
           Q(:,ii)=holder_v;
      
       end
       
   end
   end
%    En_holder=En(2);
%    En(2)=En(3);
%    En(3)=En_holder;
%    Q_holder=Q(:,2);
%    Q(:,2)=Q(:,3);
%    Q(:,3)=Q_holder;
 figure (1)  
plot3((kx)/pi,(ky)/pi,En(1)/pi,'r.','LineWidth', 1)
hold on
plot3((kx)/pi,(ky)/pi,En(2)/pi,'b.','LineWidth', 1)
hold on
plot3((kx)/pi,(ky)/pi,En(3)/pi,'g.','LineWidth', 1)
hold on
plot3((kx)/pi,(ky)/pi,En(4)/pi,'k.','LineWidth', 1)

P1=Q(:,1)*ctranspose(Q(:,1)); 
P2=Q(:,2)*ctranspose(Q(:,2)); 
P3=Q(:,3)*ctranspose(Q(:,3)); 
P4=Q(:,4)*ctranspose(Q(:,4)); 
    
    
    %%%%%%%%%%%%%%%%%%%%% dkx %%%%%%%%%%%%%
 H_1_dkx=-[0, ka.*exp(1i*(kx+dkx)),0, 0
    ka.*exp(-1i*(kx+dkx)), 0, 0, 0
     0, 0, 0, kb.*exp(1i*(kx+dkx))
     0, 0, kb.*exp(-1i*(kx+dkx)), 0];
H_2_dkx=-[0, 0,ka.*exp(1i*ky), 0
    0, 0, 0, kb.*exp(1i*ky)
     ka.*exp(-1i*ky), 0, 0, 0 
     0, kb.*exp(-1i*ky), 0, 0];
H_3_dkx=-[0, ka.*exp(-1i*(kx+dkx)),0, 0
    ka.*exp(1i*(kx+dkx)), 0, 0, 0
     0, 0, 0, kb.*exp(-1i*(kx+dkx))
     0, 0, kb.*exp(1i*(kx+dkx)), 0];
H_4_dkx=-[0, 0,ka.*exp(-1i*ky), 0
    0, 0, 0, kb.*exp(-1i*ky)
     ka.*exp(1i*ky), 0, 0, 0 
     0, kb.*exp(1i*ky), 0, 0];
 
 U1_dkx=expm(-1i.*H_1_dkx.*L./4);
 U2_dkx=expm(-1i.*H_2_dkx.*L./4);
 U3_dkx=expm(-1i.*H_3_dkx.*L./4);
 U4_dkx=expm(-1i.*H_4_dkx.*L./4);
 U_F_dkx=U4_dkx*U3_dkx*U2_dkx*U1_dkx;  
 
   
   [V_dkx,D_dkx]=eig(U_F_dkx);
   
   Q_dkx=V_dkx;

   En_dkx(1)=-log(D_dkx(1,1))./1i;
   En_dkx(2)=-log(D_dkx(2,2))./1i;
   En_dkx(3)=-log(D_dkx(3,3))./1i;
   En_dkx(4)=-log(D_dkx(4,4))./1i;
   
   for jj=1:3
   for ii=jj+1:4
       if real(En_dkx(jj))>real(En_dkx(ii))
           holder_e_dkx=En_dkx(jj);
           En_dkx(jj)=En_dkx(ii);
           En_dkx(ii)=holder_e_dkx;
           holder_v_dkx=Q_dkx(:,jj);
           Q_dkx(:,jj)=Q_dkx(:,ii);
           Q_dkx(:,ii)=holder_v_dkx;
      
       end
       
   end
   end
   
%    En_holder_dkx=En_dkx(2);
%    En_dkx(2)=En_dkx(3);
%    En_dkx(3)=En_holder_dkx;
%    Q_holder_dkx=Q_dkx(:,2);
%    Q_dkx(:,2)=Q_dkx(:,3);
%    Q_dkx(:,3)=Q_holder_dkx;
   
P1_dkx=Q_dkx(:,1)*ctranspose(Q_dkx(:,1)); 
P2_dkx=Q_dkx(:,2)*ctranspose(Q_dkx(:,2)); 
P3_dkx=Q_dkx(:,3)*ctranspose(Q_dkx(:,3)); 
P4_dkx=Q_dkx(:,4)*ctranspose(Q_dkx(:,4)); 
      

    %%%%%%%%%%%%%%%%%%%%%% dky %%%%%%%%%%%%%%%%%%%%
    
H_1_dky=-[0, ka.*exp(1i*kx),0, 0
    ka.*exp(-1i*kx), 0, 0, 0
     0, 0, 0, kb.*exp(1i*kx)
     0, 0, kb.*exp(-1i*kx), 0];
H_2_dky=-[0, 0,ka.*exp(1i*(ky+dky)), 0
    0, 0, 0, kb.*exp(1i*(ky+dky))
     ka.*exp(-1i*(ky+dky)), 0, 0, 0 
     0, kb.*exp(-1i*(ky+dky)), 0, 0];
H_3_dky=-[0, ka.*exp(-1i*kx),0, 0
    ka.*exp(1i*kx), 0, 0, 0
     0, 0, 0, kb.*exp(-1i*kx)
     0, 0, kb.*exp(1i*kx), 0];
H_4_dky=-[0, 0,ka.*exp(-1i*(ky+dky)), 0
    0, 0, 0, kb.*exp(-1i*(ky+dky))
     ka.*exp(1i*(ky+dky)), 0, 0, 0 
     0, kb.*exp(1i*(ky+dky)), 0, 0];
 
 U1_dky=expm(-1i.*H_1_dky.*L./4);
 U2_dky=expm(-1i.*H_2_dky.*L./4);
 U3_dky=expm(-1i.*H_3_dky.*L./4);
 U4_dky=expm(-1i.*H_4_dky.*L./4);
 U_F_dky=U4_dky*U3_dky*U2_dky*U1_dky;  
 
   
  [V_dky,D_dky]=eig(U_F_dky);
   
   Q_dky=V_dky;

   En_dky(1)=-log(D_dky(1,1))./1i;
   En_dky(2)=-log(D_dky(2,2))./1i;
   En_dky(3)=-log(D_dky(3,3))./1i;
   En_dky(4)=-log(D_dky(4,4))./1i;
   
   for jj=1:3
   for ii=jj+1:4
       if real(En_dky(jj))>real(En_dky(ii))
           holder_e_dky=En_dky(jj);
           En_dky(jj)=En_dky(ii);
           En_dky(ii)=holder_e_dky;
           holder_v_dky=Q_dky(:,jj);
           Q_dky(:,jj)=Q_dky(:,ii);
           Q_dky(:,ii)=holder_v_dky;
      
       end
       
   end
   end
   
%    En_holder_dky=En_dky(2);
%    En_dky(2)=En_dky(3);
%    En_dky(3)=En_holder_dky;
%    Q_holder_dky=Q_dky(:,2);
%    Q_dky(:,2)=Q_dky(:,3);
%    Q_dky(:,3)=Q_holder_dky;
   
P1_dky=Q_dky(:,1)*ctranspose(Q_dky(:,1)); 
P2_dky=Q_dky(:,2)*ctranspose(Q_dky(:,2)); 
P3_dky=Q_dky(:,3)*ctranspose(Q_dky(:,3)); 
P4_dky=Q_dky(:,4)*ctranspose(Q_dky(:,4));    

    %%%%%%%%%%%%%%%%
    A1kx=(P1_dkx-P1)./dkx;
    A2kx=(P2_dkx-P2)./dkx;
    A3kx=(P3_dkx-P3)./dkx;
    A4kx=(P4_dkx-P4)./dkx;
    
    A1ky=(P1_dky-P1)./dky;
    A2ky=(P2_dky-P2)./dky;
    A3ky=(P3_dky-P3)./dky;
    A4ky=(P4_dky-P4)./dky;
    
    
    
    B1=A1kx*A1ky-A1ky*A1kx;
     B2=A2kx*A2ky-A2ky*A2kx;
      B3=A3kx*A3ky-A3ky*A3kx;
       B4=A4kx*A4ky-A4ky*A4kx;
   
    
   
    A1=P1*B1;
    A2=P2*B2;
    A3=P3*B3;
    A4=P4*B4;
    
    
    At1(jjx,jjy)=trace(A1);
    At2(jjx,jjy)=trace(A2);
    At3(jjx,jjy)=trace(A3);
    At4(jjx,jjy)=trace(A4);

  
     jjx=jjx+1;
    end
        jjy=jjy+1;
end


kx=(k_start*pi:stepkx:k_end*pi);
ky=(k_start*pi:stepky:k_end*pi);


I1= (trapz(ky,trapz(kx,At1,2)))/(2i*pi)
I2= (trapz(ky,trapz(kx,At2,2)))/(2i*pi)
I3= (trapz(ky,trapz(kx,At3,2)))/(2i*pi)
I4= (trapz(ky,trapz(kx,At4,2)))/(2i*pi)


%  I= (trapz(kx,At));
%  W= (trapz(ky,trapz(kx,I,2)))/(8*(pi^2))


