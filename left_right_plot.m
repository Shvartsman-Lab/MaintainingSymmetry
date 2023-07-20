clear all;
foldername='005';
load (sprintf('%s/velocityField_2D_%s.mat',foldername,foldername));
load(sprintf('g1_%s.mat',foldername));
load(sprintf('%s_midline.mat',foldername));
load(sprintf('%s/midline.mat',foldername));
vid=VideoWriter(sprintf('%s/left_right_%s.mp4',foldername,foldername),'MPEG-4');
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
vid.FrameRate=7;
open(vid);
time=size(u_original,1);
load(sprintf('%s/analysis.mat',foldername));
for (i=1:181)
    z_dash=para1.zmin+(para1.zmax-para1.zmin)*(i-1)/180;
    R_p=sqrt(polyval (para1.pR,z_dash,para1.SR,para1.muR));
    e_p=polyval (para1.pe,z_dash,para1.Se,para1.mue);
    R_p_2=R_p/sqrt (1-e_p*e_p);
    x_0=polyval (para1.pX,z_dash,para1.SX,para1.muX);
    y_0=polyval (para1.pY,z_dash,para1.SY,para1.muY);
    for (j=1:81)
        phi_h=-pi+(j-1)*pi/40;
        Z(i,j)=z_dash;
        X(i,j)=R_p*cos (phi_h+para1.pphase)+x_0;
        Y(i,j)=R_p_2*sin(phi_h+para1.pphase)+y_0;
    end
end
[Nx,Ny,Nz] = surfnorm(X,Y,Z);
Z_a=para1.zmin+(para1.zmax-para1.zmin)*(0:1/180:1);
phi_a=-pi:pi/40:pi;

vid=VideoWriter(sprintf('%s/L_R_V_%s.mp4',foldername,foldername),'MPEG-4');
vid.FrameRate=7;
open(vid);
for (t=6:time-t_start-1)
    surf (0.98*(Z-shift*Nz),0.98*(Y-shift*Ny),0.98*(X-shift*Nx),'edgecolor','none','Facecolor',[0.7 0.7 0.7]);
    axis equal;
    axis off;
    hold on; 
    phi_midline=nanmean(phi(1:300,t_start+t));
    for (i=1:size(z_range,2)-10)
        z_midline(i)=z_range(i);
        z_ots=interp1(Z(:,41)-shift*Nz(:,41),Z(:,41),z_range(i));
        R_p=sqrt(polyval (para1.pR,z_ots,para1.SR,para1.muR));
        e_p=polyval (para1.pe,z_ots,para1.Se,para1.mue);
        R_p_2=R_p/sqrt(1-e_p*e_p);
        x_0=polyval (para1.pX,z_ots,para1.SX,para1.muX);
        y_0=polyval (para1.pY,z_ots,para1.SY,para1.muY);
        Nx_1=interp2(phi_a,Z_a,Nx,-para1.pphase,double(z_ots));
        Ny_2=interp2(phi_a,Z_a,Ny,pi/2-para1.pphase,double(z_ots));
        R_p=R_p-shift*Nx_1;
        R_p_2=R_p_2-shift*Ny_2;
        x_midline(i)=x_0+R_p*cos(phi_midline+para1.pphase);
        y_midline(i)=y_0+R_p_2*sin(phi_midline+para1.pphase);
    end
    load(sprintf('%s/analysis.mat',foldername));
    z_grid=linspace(min(coord(:,1)),max(coord(:,1)),20);
    phi_grid=-pi:pi/20:pi;
    [X_2,Y_2]=meshgrid(phi_grid,double(z_grid));
    for (i=1:size(coord,1))
        coord(i,2)=coord(i,2)-phi_midline;
        if (coord(i,2)>-pi) coord(i,2)=coord(i,2)+2*pi ; end
        if (coord(i,2)>pi) coord(i,2)=coord(i,2)-2*pi; end
    end
    l=size(coord,1);
    extra=l;
    for (i=1:l)
        if (coord(i,2)>3*pi/4)
            extra=extra+1;
            coord(extra,1)=coord(i,1); 
            coord(extra,2)=coord(i,2)-2*pi;
            v_theta(extra,:)=v_theta(i,:);
            v_phi(extra,:)=v_phi(i,:);
        end
        if (coord(i,2)<-3*pi/4)
            extra=extra+1;
            coord(extra,1)=coord(i,1); 
            coord(extra,2)=coord(i,2)+2*pi;
            v_theta(extra,:)=v_theta(i,:);
            v_phi(extra,:)=v_phi(i,:);
        end
    end
    F_theta=scatteredInterpolant(double(coord(:,1)),double(coord(:,2)),double(v_theta(:,t_start+t)),'nearest','none');
    F_phi=scatteredInterpolant(double(coord(:,1)),double(coord(:,2)),double(v_phi(:,t_start+t)),'nearest','none');
    v_theta_grid=F_theta(Y_2,X_2);
    v_phi_grid=F_phi(Y_2,X_2);
    v_phi_grid_image=F_phi(Y_2,-X_2);
    for (j=1:size(X_2,1))
        z_h=Y_2(j,1);
        if (isnan (z_act(j))) continue; end ;
        z_ots=interp1(Z(:,41)-shift*Nz(:,41),Z(:,41),z_h);
        R_p(j)=sqrt(polyval (para1.pR,z_ots,para1.SR,para1.muR));
        e_p=polyval (para1.pe,z_ots,para1.Se,para1.mue);
        R_p_2(j)=R_p(j)/sqrt(1-e_p*e_p);
        x_0(j)=polyval (para1.pX,z_ots,para1.SX,para1.muX);
        y_0(j)=polyval (para1.pY,z_ots,para1.SY,para1.muY);
        Nx_1=interp2(phi_a,Z_a,Nx,-para1.pphase,double(z_ots));
        Ny_2=interp2(phi_a,Z_a,Ny,pi/2-para1.pphase,double(z_ots));
        R_p(j)=R_p(j)-shift*Nx_1;
        R_p_2(j)=R_p_2(j)-shift*Ny_2;
        d_x_0_d_z(j)=polyval (polyder(para1.pX),z_ots,[],para1.muX)/para1.muX(2);
        d_y_0_d_z(j)=polyval (polyder(para1.pY),z_ots,[],para1.muY)/para1.muY(2);
    end 
    
    for (i=2:size(X_2,1)-1)
        z_h=Y_2(i,1);
        d_R_p_d_z(i)=(R_p(i+1)-R_p(i-1))/(Y_2(i+1,1)-Y_2(i-1,1));
        d_R_p_2_d_z(i)=(R_p_2(i+1)-R_p_2(i-1))/(Y_2(i+1,1)-Y_2(i-1,1));
        f_z(:,1)=fp(:,floor (z_h-z_range(1)+1))*(1-(z_h-floor (z_h)))+fp(:,floor (z_h)-z_range(1)+2)*(z_h-floor(z_h));
        gff_int(:,1)=gff(:,floor (z_h-z_range(1)+1))*(1-(z_h-floor (z_h)))+fp(:,floor (z_h)-z_range(1)+2)*(z_h-floor (z_h));
        for (j=1:size(X_2,2))       
            phi_h=X_2(i,j);
            Z_3d(i,j)=z_h;
            R_p_h=R_p(i); R_p_2_h=R_p_2(i);
            X_3d(i,j)=x_0(i)+R_p_h*cos (phi_h+phi_midline+para1.pphase);
            Y_3d(i,j)=y_0(i)+R_p_2_h*sin(phi_h+phi_midline+para1.pphase);
            d_q_d_z=interp1(z_range,sqrt(gzz(II,:)),z_h);
            V_Z(i,j)=v_theta_grid(i,j)/d_q_d_z;
            d_phi_d_t= v_phi_grid(i,j)/sqrt((R_p_h*sin(phi_h+phi_midline+para1.pphase))^2+(R_p_2_h*cos(phi_h+phi_midline+para1.pphase))^2);
            V_X(i,j)=-R_p_h*sin(phi_h+phi_midline+para1.pphase)*d_phi_d_t+(d_x_0_d_z(i)+d_R_p_d_z(i)*cos(phi_h+phi_midline+para1.pphase))*v_z(i,j);
            V_Y(i,j)=R_p_2_h*cos(phi_h+phi_midline+para1.pphase)*d_phi_d_t+(d_y_0_d_z(i)+d_R_p_2_d_z(i)*sin(phi_h+phi_midline+para1.pphase))*v_z(i,j);  
            d_phi_d_t=v_phi_grid_image(i,j)/sqrt((R_p_h*sin(phi_h+phi_midline+para1.pphase))^2+(R_p_2_h*cos(phi_h+phi_midline+para1.pphase))^2);
            V_X_2(i,j)=-R_p_h*sin(phi_h+phi_midline+para1.pphase)*d_phi_d_t+(d_x_0_d_z(i)+d_R_p_d_z(i)*cos(phi_h+phi_midline+para1.pphase))*v_z(i,j);
            V_Y_2(i,j)=R_p_2_h*cos(phi_h+phi_midline+para1.pphase)*d_phi_d_t+(d_y_0_d_z(i)+d_R_p_2_d_z(i)*sin(phi_h+phi_midline+para1.pphase))*v_z(i,j);                  
        end
    end
    %scatter3(Z_3d(:),Y_3d(:),X_3d(:),4,'b','filled');
    plot3(z_midline,y_midline,x_midline,'color','k','Linewidth',1.2);
    quiver3(Z_3d,Y_3d,X_3d,2*V_Z,2*V_Y,2*V_X,0,'color','b','Linewidth',0.3,'MaxHeadSize',50);
     quiver3(Z_3d,Y_3d,X_3d,2*V_Z,-2*V_Y_2,-2*V_X_2,0,'color','r','Linewidth',0.3,'MaxHeadSize',50);
    view ([0 90+(phi_midline+para1.pphase)*180/pi]);
    F=getframe(gcf);
    writeVideo(vid,F);
    hold off;
    
end
 close (vid);



% 
% 
% 
% 
% 
% %%
% 
% phi_interp(isnan(phi_interp))=[];
% phi_midline=mean(phi_interp);
% 
% for (i=1:size(coord,1))
%     coord(i,2)=coord(i,2)-phi_midline;
%     if (coord(i,2)>-pi) coord(i,2)=coord(i,2)+2*pi ; end
%     if (coord(i,2)>pi) coord(i,2)=coord(i,2)-2*pi; end
% end
% 
% %%
% pos_p=y{1};
% pos_q=x{1};
% sz=size(pos_p);
% z_range=double (para1.zmin:dz:para1.zmax);
% mid=size(z_range,2);shift=(size(z_range,2)-size(zp,2))/2+1;
% z_range(1:shift)=[];z_range(size(zp,2)+1:end)=[];
% df=double(df);
% phi_range=-pi:df:pi;
% fp=fp+(bdry{2}(2)-bdry{2}(1))/2;
% id=0;
% d_p_d_z_set=zeros (size(zp));
% zp = real(cumsum(sqrt(gzz)*dz, 2));
% [~,maxind] =  max(zp(:));
% [II,~] = ind2sub(size(zp),maxind);
% d_p_d_z_set(:,1:end-1)=(fp(:,2:end)-fp(:,1:end-1))/dz;
% 
% 
% 
% for (t=1:time)
%      
%      v_phi_grid_2=F_phi(Y,-X);
%      quiver(q,p,0.5*v_theta_grid,0.5*v_phi_grid,0);
%      hold on
%      quiver (q,p,0.5*v_theta_grid,-0.5*v_phi_grid_2,0);
%      plot (q(:,1),mid_phi,'Linestyle','--');
%      axis equal;
%      xlim ([0 1000]);
%      ylim ([0 1000]);
%      F=getframe(gcf);
%      writeVideo(vid,F);
%      hold off;
% end
% close(vid);
% 
% 
% 
