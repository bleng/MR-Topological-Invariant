close all
clear 
clc
%%%%%%%%%%%%%%%
L=1;
thetaa=0.45*pi;
thetab=0.01*pi;

Bg1=-0.2*pi/L;
Bg3=-Bg1;
Bg2=1*pi/L;

Bg_test=Bg2;

ka=thetaa/(L/4);
kb=thetab/(L/4);

stepz=0.1*L;
stepkx=0.01*pi;
stepky=0.01*pi;
dz=0.001*stepz;
dkx=0.001*stepkx;
dky=0.001*stepky;


k_start=-0.5;
k_end=0.5;

kx=(k_start*pi:stepkx:k_end*pi);
ky=(k_start*pi:stepky:k_end*pi);
Z=0:stepz:L;



At=zeros(length(kx),length(ky), length(Z));

jjx=1;
jjy=1;
jjz=1;
for Z=0:stepz:L
    %dz=0.1*stepz;
    jjy=1;
for ky=(k_start*pi:stepkx:k_end*pi)
    jjx=1;
    for kx=(k_start*pi:stepky:k_end*pi)
        
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
 
 %%%%%%%%%%% U(k,z) %%%%%%%%%%%
  
     
     if Z<=(L/4)
       U_z=expm(-1i.*H_1.*Z);
     elseif (Z>(L/4)) && (Z<=(L/2))
          U_z=expm(-1i.*H_2.*(Z-(L/4)))*U1;
     elseif (Z>(L/2)) && (Z<=(3*L/4))
         U_z=expm(-1i.*H_3.*(Z-(L/2)))*U2*U1;
     elseif (Z>(3*L/4)) && (Z<=L)
         U_z=expm(-1i.*H_4.*(Z-(3*L/4)))*U3*U2*U1;
     end
    
    
 

   
   
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
   
   
   H_eff_e=diag(1*En./L);
   
   for ii=1:4
       if H_eff_e(ii,ii)<Bg_test
           H_eff_e(ii,ii)=H_eff_e(ii,ii)+(2*pi/L);
       end
   end

    V_eps=Q*(expm(1i.*H_eff_e.*Z))*inv(Q);
   
    
   
    U_eps=U_z*V_eps;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%% dz %%%%%%%%%%5
    
   % U_z_d=zeros(4,4);
    if Z+dz<=L/4
       U_z_d=expm(-1i.*H_1.*(Z+dz));
     elseif (Z+dz>(L/4)) && (Z+dz<=(L/2))
          U_z_d=expm(-1i.*H_2.*(Z+dz-(L/4)))*U1;
     elseif (Z+dz>(L/2)) && (Z+dz<=(3*L/4))
         U_z_d=expm(-1i.*H_3.*(Z+dz-(L/2)))*U2*U1;
     elseif (Z+dz>(3*L/4)) && (Z+dz<=L)
         U_z_d=expm(-1i.*H_4.*(Z+dz-(3*L/4)))*U3*U2*U1;
         elseif Z+dz>L
              U_z_d=expm(-1i.*H_4.*(Z-dz-(3*L/4)))*U3*U2*U1;
      %        U_z_d=expm(-1i.*H_1.*(Z+dz-L));
              dz=-dz;
             
             
    end
     
  %  V_eps_dz=Q*(expm(-ang_U.*(Z+dz)))*inv(Q);
   %  V_eps_dz=expm(1i*H_eff*(Z+dz));
    V_eps_dz=Q*(expm(1i.*H_eff_e.*(Z+dz)))*inv(Q);
    
    U_eps_dz=U_z_d*V_eps_dz;
   
    du_eps_z=(U_eps_dz-U_eps)./dz;
    
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
   
   
   H_eff_e_dkx=diag(1*En_dkx./L);
   
   for ii=1:4
       if H_eff_e_dkx(ii,ii)<Bg_test
           H_eff_e_dkx(ii,ii)=H_eff_e_dkx(ii,ii)+(2*pi/L);
       end
   end

    V_eps_dkx=Q_dkx*(expm(1i.*H_eff_e_dkx.*Z))*inv(Q_dkx);
  
   

   
%     V_eps_dkx=Q_dkx*(expm(1i*ang_U_dkx.*Z))*inv(Q_dkx);
    
    
    % U_z_dkx=zeros(4,4);
     
     
     if Z<=L/4
       U_z_dkx=expm(-1i.*H_1_dkx.*Z);
     elseif (Z>(L/4)) && (Z<=(L/2))
          U_z_dkx=expm(-1i.*H_2_dkx.*(Z-(L/4)))*U1_dkx;
     elseif (Z>(L/2)) && (Z<=(3*L/4))
         U_z_dkx=expm(-1i.*H_3_dkx.*(Z-(L/2)))*U2_dkx*U1_dkx;
     elseif (Z>(3*L/4)) && (Z<=L)
         U_z_dkx=expm(-1i.*H_4_dkx.*(Z-(3*L/4)))*U3_dkx*U2_dkx*U1_dkx;
     end
    
    
    U_eps_dkx=U_z_dkx*V_eps_dkx;
    
    du_eps_kx=(U_eps_dkx-U_eps)./dkx;
    
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
   
   
   H_eff_e_dky=diag(1*En_dky./L);
   for ii=1:4
       if H_eff_e_dky(ii,ii)<Bg_test
           H_eff_e_dky(ii,ii)=H_eff_e_dky(ii,ii)+(2*pi/L);
       end
   end

    V_eps_dky=Q_dky*(expm(1i.*H_eff_e_dky.*Z))*inv(Q_dky);
   
  
   % V_eps_dky=Q_dky*(expm(1i*ang_U_dky.*Z))*inv(Q_dky);
    
    
     %U_z_dky=zeros(4,4);
     
     
     if Z<=L/4
       U_z_dky=expm(-1i.*H_1_dky.*Z);
     elseif (Z>(L/4)) && (Z<=(L/2))
          U_z_dky=expm(-1i.*H_2_dky.*(Z-(L/4)))*U1_dky;
     elseif (Z>(L/2)) && (Z<=(3*L/4))
         U_z_dky=expm(-1i.*H_3_dky.*(Z-(L/2)))*U2_dky*U1_dky;
     elseif (Z>(3*L/4)) && (Z<=L)
         U_z_dky=expm(-1i.*H_4_dky.*(Z-(3*L/4)))*U3_dky*U2_dky*U1_dky;
     end
    
    
    U_eps_dky=U_z_dky*V_eps_dky;
    
    du_eps_ky=(U_eps_dky-U_eps)./dky;
    %%%%%%%%%%%%%%%%
    
    Akx=inv(U_eps)*du_eps_kx;
    Aky=inv(U_eps)*du_eps_ky;
    Az=inv(U_eps)*du_eps_z;
    
    B=Akx*Aky-Aky*Akx;
    A=Az*B;
    At(jjx,jjy,jjz)=trace(A);
%     if abs(trace(B))>1e-5
%         trace(B)
%         jjz
%     end
    
   %Sum=Sum+At(jjx,jjy,jjz); 
    
    
    
    
    
    
    
   

        
   
  
     jjx=jjx+1;
    end
        jjy=jjy+1;
end
    jjz=jjz+1;
end

kx=(k_start*pi:stepkx:k_end*pi);
ky=(k_start*pi:stepky:k_end*pi);
Z=0:stepz:L;

I= (trapz(ky,trapz(kx,At,2)));
W= (trapz(Z,I))/(8*(pi^2))

%  I= (trapz(kx,At));
%  W= (trapz(ky,trapz(kx,I,2)))/(8*(pi^2))


%W=Sum*stepz*stepkx*stepky/(8*pi^2)
%W=Sum*stepz*stepkx*stepky/(8*pi^2)