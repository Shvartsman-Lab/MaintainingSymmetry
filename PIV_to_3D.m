clear all;
cd WT;
foldername='092639';
t_start=6;
load (sprintf('%s/velocityField_2D_%s.mat',foldername,foldername));
load(sprintf('%s/g1_%s.mat',foldername,foldername));
% u_original=u_filtered;
% v_original=v_filtered;
time=size(u_original,1);
pos_p=y{1};
pos_q=x{1};
sz=size(pos_p);
z_range=double (para1.zmin:dz:para1.zmax);
mid=size(z_range,2);shift=(size(z_range,2)-size(zp,2))/2+1;
z_range(1:shift)=[];z_range(size(zp,2)+1:end)=[];
df=double(df);
phi_range=-pi:df:pi;
fp=fp+(bdry{2}(2)-bdry{2}(1))/2;
id=0;
d_p_d_z_set=zeros (size(zp));
zp = real(cumsum(sqrt(gzz)*dz, 2));
[~,maxind] =  max(zp(:));
[II,~] = ind2sub(size(zp),maxind);
d_p_d_z_set(:,1:end-1)=(fp(:,2:end)-fp(:,1:end-1))/dz;
% im=imread (sprintf('%s/cylinder1_proper_%s/cmp_1_1_T0009.tif',foldername,foldername));
% imshow (im,[0 max(max(im))])
% hold on;
% plot (zp(1,:),fp (floor (size(fp,1)),:),'Linewidth',2);
% hold off
for (i=1:181)
    z_dash=para1.zmin+(para1.zmax-para1.zmin)*(i-1)/180;
    R_p=sqrt(polyval (para1.pR,z_dash,para1.SR,para1.muR));
    e_p=polyval (para1.pe,z_dash,para1.Se,para1.mue);
    R_p_2=R_p/sqrt (1-e_p*e_p);
    x_0=polyval (para1.pX,z_dash,para1.SX,para1.muX);
    y_0=polyval (para1.pY,z_dash,para1.SY,para1.muY);
    for (j=1:81)
        phi=-pi+(j-1)*pi/40;
        Z(i,j)=z_dash;
        X(i,j)=R_p*cos (phi+para1.pphase)+x_0;
        Y(i,j)=R_p_2*sin(phi+para1.pphase)+y_0;
    end
end
[Nx,Ny,Nz] = surfnorm(X,Y,Z);
Z_a=para1.zmin+(para1.zmax-para1.zmin)*(0:1/180:1);
phi_a=-pi:pi/40:pi;

%%

for (j=1:sz(2))
    z=pos_q(1,j);
    z_act(j)=interp1(zp(1,:),z_range,z-bdry{1}(1),'spline',0/0);
    if (isnan (z_act(j))) continue; end ;
    z_ots=interp1(Z(:,41)-shift*Nz(:,41),Z(:,41),z_act(j));
    R_p(j)=sqrt(polyval (para1.pR,z_ots,para1.SR,para1.muR));
    e_p=polyval (para1.pe,z_ots,para1.Se,para1.mue);
    R_p_2(j)=R_p(j)/sqrt(1-e_p*e_p);
    x_0(j)=polyval (para1.pX,z_ots,para1.SX,para1.muX);
    y_0(j)=polyval (para1.pY,z_ots,para1.SY,para1.muY);
    Nx_1=interp2(phi_a,Z_a,Nx,-para1.pphase,double(z_ots));
    Ny_2=interp2(phi_a,Z_a,Ny,pi/2-para1.pphase,double(z_ots));
    R_p(j)=R_p(j)-shift*Nx_1;
    R_p_2(j)=R_p_2(j)-shift*Ny_2;
    d_x_0_d_z(j)=polyval (polyder(para1.pX),z_act(j),[],para1.muX)/para1.muX(2);
    d_y_0_d_z(j)=polyval (polyder(para1.pY),z_act(j),[],para1.muY)/para1.muY(2);
end
id=0;
for(j=2:sz(2)-1)
    d_R_p_d_z=(R_p(j+1)-R_p(j-1))/(z_act(j+1)-z_act(j-1));
    d_R_p_2_d_z=(R_p_2(j+1)-R_p_2(j-1))/(z_act(j+1)-z_act(j-1));
    f_z(:,1)=fp(:,floor (z_act(j)-z_range(1)+1))*(1-(z_act(j)-floor (z_act(j))))+fp(:,floor (z_act(j))-z_range(1)+2)*(z_act(j)-floor (z_act(j)));
    gff_int(:,1)=gff(:,floor (z_act(j)-z_range(1)+1))*(1-(z_act(j)-floor (z_act(j))))+fp(:,floor (z_act(j))-z_range(1)+2)*(z_act(j)-floor (z_act(j)));
    for (i=1:sz(1))
            phi =interp1(f_z,phi_range,pos_p(i,j),'spline',0/0);
            if (isnan (phi)) continue; end ;
            id=id+1;
            coord(id,1)=z_act(j);
            coord(id,2)=phi;
            coord_3d(id,3)=z_act(j);
            coord_3d(id,1)=x_0(j)+R_p(j)*cos (coord(id,2)+para1.pphase);
            coord_3d(id,2)=y_0(j)+R_p_2(j)*sin(coord(id,2)+para1.pphase);
            for (t=1:time)
                v_p=v_original{t};
                v_q=u_original{t};
                if (isnan(v_p(i,j))) v_p(i,j)=0; end; 
                if (isnan(v_q(i,j))) v_q(i,j)=0; end; 
                v_theta(id,t)=v_q(i,j);
                d_q_d_z=interp1(z_range,sqrt(gzz(II,:)),z_act(j));
                v_z(id,t)=v_q(i,j)/d_q_d_z;
                d_p_d_phi=interp1(phi_range,sqrt(gff_int),phi);
                d_p_d_z=interp1(phi_range,d_p_d_z_set(:,floor(z_act(j)-z_range(1)+1)),phi);
                d_phi_d_z=-d_p_d_z/d_p_d_phi;
                d_phi_d_t= d_phi_d_z*v_z(id,t)+v_p(i,j)/d_p_d_phi;
                d_R_p_d_z_h(id)=d_R_p_d_z;
                d_R_p_2_d_z_h(id)=d_R_p_2_d_z;
                v_phi(id,t)=sqrt ((R_p(j)*sin(phi+para1.pphase))^2+(R_p_2(j)*cos(phi+para1.pphase))^2)*d_phi_d_t;
                v_x(id,t)=-R_p(j)*sin (coord(id,2)+para1.pphase)*d_phi_d_t+(d_x_0_d_z(j)+d_R_p_d_z*cos (coord(id,2)+para1.pphase))*v_z(id,t);
                v_y(id,t)=R_p_2(j)*cos (coord(id,2)+para1.pphase)*d_phi_d_t+(d_y_0_d_z(j)+d_R_p_2_d_z*sin(coord(id,2)+para1.pphase))*v_z(id,t);
    end
        disp (j);
    end
end

 
%% plotting
vid=VideoWriter(sprintf('%s/timeseries_%s.mp4',foldername,foldername),'MPEG-4');
vid.FrameRate=7;
open(vid);
for (t=1:time)
    surf (Z-shift*Nz,Y-shift*Ny,X-shift*Nx,'edgecolor','none','Facecolor',[0.7 0.7 0.7]);
    axis equal;
    axis off;
    hold on;
    quiver3(coord_3d(1:id,3),coord_3d(1:id,2),coord_3d(1:id,1),v_z(1:id,t),v_y(1:id,t),v_x(1:id,t),'Color','b','Linewidth',1.5,'Maxheadsize',0.1);
    %scatter3(coord_3d(1:id,1),coord_3d(1:id,2),coord_3d(1:id,3),40,'r','filled');
    view ([10 0]);
    %hold off;
    F=getframe(gcf);
    writeVideo(vid,F);
    hold off;
end
 close(vid);
%% SVD

A=cat(1,v_theta (1:id,t_start:t_start+60),v_phi (1:id,t_start:t_start+60)); 
[U,S,V]=svd(A-mean(A,1),'econ'); 
E=S*S;
%U(:,1)=-U(:,1);
V(:,1)=-V(:,1);
%%
for (mode=1:3)
    figure (1);
    %subplot(2,3,mode);%, 'sh', 0.04, 'sv', 0.05, 'padding', 0, 'margin', 0);
    plot (0:60,V(:,mode),'k','LineWidth',1.2);
    x=xlabel ('$t$ (min)','Interpreter','latex','FontSize',25);
    xticks (0:20:100);
    set(gca,'linewidth',0.2);
    set(gca,'ticklabelinterpreter','latex');
    set(gca,'Fontsize',25);
    axis ([ 0 98 -0.5 0.5]);
    box on;
    hold off
    saveas (gcf,strcat (foldername,'/time_mode_',num2str(mode),'.fig'));
    
    figure (2);
    surf (0.98*(Z-shift*Nz),0.98*(Y-shift*Ny),0.98*(X-shift*Nx),'edgecolor','none','Facecolor',[0 0 0]);
    axis equal;
    axis off;
    hold on;
    v_mode=zeros (id,3);
    
    for (id=1:size(coord,1))
        z_h=coord(id,1); 
        phi=coord(id,2);
        z_ots=interp1(Z(:,41)-shift*Nz(:,41),Z(:,41),z_h);
        R_p_h=sqrt(polyval (para1.pR,z_ots,para1.SR,para1.muR));
        e_p=polyval (para1.pe,z_ots,para1.Se,para1.mue);
        R_p_2_h=R_p_h/sqrt(1-e_p*e_p);
        x_0=polyval (para1.pX,z_ots,para1.SX,para1.muX);
        y_0=polyval (para1.pY,z_ots,para1.SY,para1.muY);
        Nx_1=interp2(phi_a,Z_a,Nx,-para1.pphase,double(z_ots));
        Ny_2=interp2(phi_a,Z_a,Ny,pi/2-para1.pphase,double(z_ots));
        R_p_h=R_p_h-shift*Nx_1;
        R_p_2_h=R_p_2_h-shift*Ny_2;
        v_mode(id,3)=U(id,mode);
        d_q_d_z=interp1(z_range,sqrt(gzz(II,:)),z_h);
        v_mode(id,3)=v_mode(id,3)/d_q_d_z;
        v_phi=U(545+id,mode);
        d_phi_d_t= v_phi/sqrt((R_p_h*sin(phi+para1.pphase))^2+(R_p_2_h*cos(phi+para1.pphase))^2);
        d_x_0_d_z_h=polyval(polyder(para1.pX),z_h,[],para1.muX)/para1.muX(2);
        d_y_0_d_z_h=polyval(polyder(para1.pY),z_h,[],para1.muY)/para1.muY(2);
        v_mode(id,1)=(d_x_0_d_z_h+d_R_p_d_z_h(id)*cos(phi+para1.pphase))*v_mode(id,3)-R_p_h*sin(phi+para1.pphase)*d_phi_d_t;
        v_mode(id,2)=(d_y_0_d_z_h+d_R_p_2_d_z_h(id)*sin(phi+para1.pphase))*v_mode(id,3)+R_p_2_h*cos(phi+para1.pphase)*d_phi_d_t;               
        
    end
    quiver3(coord_3d(1:id,3),coord_3d(1:id,2),coord_3d(1:id,1),v_mode(1:id,3),v_mode(1:id,2),v_mode(1:id,1),'Color',[0.75 0.75 0],'Linewidth',1.5,'Maxheadsize',0.5);  
    view ([180 0]); hold off;
    saveas (gcf,strcat (foldername,'/space_mode_',num2str(mode),'.fig'));
    hold off;
    disp (mode);
  
end
save (sprintf ('%s/analysis.mat',foldername),'coord','coord_3d','v_x','v_y','v_z','v_theta','v_phi','U','S','V');

