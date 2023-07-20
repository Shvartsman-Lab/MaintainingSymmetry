clear all;
foldername='005';
t_start=6;
load (sprintf('%s/velocityField_2D_%s.mat',foldername,foldername));
load(sprintf('g1_%s.mat',foldername));
load (sprintf ('%s/analysis.mat',foldername));
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


mode=1;
l=size(coord,1);
extra=l;
for (i=1:l)
    d_q_d_z=interp1(z_range,sqrt(gzz(II,:)),coord(i,1));
    v_z_mode(i,1)=U(i,mode)/d_q_d_z;
    v_phi_mode(i,1)=U(i+l,mode);
    if (coord(i,2)>3*pi/4)
            extra=extra+1;
            coord(extra,1)=coord(i,1); 
            coord(extra,2)=coord(i,2)-2*pi;
            v_z_mode(extra,:)=v_z_mode(i,:);
            v_phi_mode(extra,:)=v_phi_mode(i,:);
    end
     if (coord(i,2)<-3*pi/4)
            extra=extra+1;
            coord(extra,1)=coord(i,1); 
            coord(extra,2)=coord(i,2)+2*pi;
            v_z_mode(extra,:)=v_z_mode(i,:);
            v_phi_mode(extra,:)=v_phi_mode(i,:);
     end
end
a=[ones(13,1)*min(z_range) (-1.2*pi:pi/5:1.2*pi)'];
b=[ones(13,1)*max(z_range) (-1.2*pi:pi/5:1.2*pi)'];
coord=[coord; a;b];
v_z_mode=[v_z_mode; zeros(26,1)];
v_phi_mode=[v_phi_mode;zeros(26,1)];
F_z=scatteredInterpolant(double(coord(:,1)),double(coord(:,2)),double(v_z_mode),'natural','none');
F_phi=scatteredInterpolant(double(coord(:,1)),double(coord(:,2)),double(v_phi_mode),'natural','none');
%%
n_seeds=200;
t_max=10000;
stream_coord=zeros (n_seeds,t_max,2);
stream_coord_3d=zeros (n_seeds,t_max,3);
% for (i=1:n_seeds)
%     stream_coord(i,1,1)=min(coord(:,1))+rand()*(max(coord(:,1))-min(coord(:,1)));
%     stream_coord(i,1,2)=-pi+rand()*2*pi;
% end
[Z_g,phi_g]=meshgrid(linspace(0.9*min(coord(:,1)),0.9*max(coord(:,1)),20),(-0.8*pi:pi/5:pi));
stream_coord(:,1,1)=Z_g(:);
stream_coord(:,1,2)=phi_g(:);
flag=0;
for (t=1:t_max-1)
    for (i=1:n_seeds)
        z=stream_coord(i,t,1);
        phi=stream_coord(i,t,2);
        z_ots=interp1(Z(:,41)-shift*Nz(:,41),Z(:,41),z);
        R_p_h=sqrt(polyval (para1.pR,z_ots,para1.SR,para1.muR));
        e_p=polyval (para1.pe,z_ots,para1.Se,para1.mue);
        R_p_2_h=R_p_h/sqrt(1-e_p*e_p);
        x_0=polyval (para1.pX,z_ots,para1.SX,para1.muX);
        y_0=polyval (para1.pY,z_ots,para1.SY,para1.muY);
        Nx_1=interp2(phi_a,Z_a,Nx,-para1.pphase,double(z_ots));
        Ny_2=interp2(phi_a,Z_a,Ny,pi/2-para1.pphase,double(z_ots));
        R_p_h=R_p_h-shift*Nx_1;
        R_p_2_h=R_p_2_h-shift*Ny_2;
        stream_coord_3d(i,t,3)=z;
        stream_coord_3d(i,t,1)= x_0+R_p_h*cos(phi+para1.pphase);
        stream_coord_3d(i,t,2)= y_0+R_p_2_h*sin(phi+para1.pphase);
        if (isnan(stream_coord_3d(i,t,1))) disp(i); flag=1; break; end
        v_1=F_z(z,phi);
        v_2=F_phi(z,phi)/sqrt ((R_p_h*sin(phi+para1.pphase))^2+(R_p_2_h*cos(phi+para1.pphase))^2);
        if (isnan(v_2)) disp(i); flag=1; break; end
        stream_coord(i,t+1,1)=stream_coord(i,t,1)+v_1;
        stream_coord(i,t+1,2)=stream_coord(i,t,2)+v_2;
        if (stream_coord(i,t+1,2)>pi) stream_coord(i,t+1,2)=stream_coord(i,t+1,2)-2*pi; end
        if (stream_coord(i,t+1,2)<-pi) stream_coord(i,t+1,2)=stream_coord(i,t+1,2)+2*pi; end
        if (stream_coord(i,t+1,1)<min(0.99*z_range))
            stream_coord(i,t+1,1)=2*min(0.99*z_range)-stream_coord(i,t+1,1);
            stream_coord(i,t+1,2)=stream_coord(i,t+1,2)+pi;
            if (stream_coord(i,t+1,2)>pi) stream_coord(i,t+1,2)=stream_coord(i,t+1,2)-2*pi; end
        end
        if (stream_coord(i,t+1,1)>max(0.99*z_range))
            stream_coord(i,t+1,1)=2*max(0.99*z_range)-stream_coord(i,t+1,1);
            stream_coord(i,t+1,2)=stream_coord(i,t+1,2)+pi;
            if (stream_coord(i,t+1,2)>pi) stream_coord(i,t+1,2)=stream_coord(i,t+1,2)-2*pi; end
        end
        
    end
    disp(t);
    if (flag) break; end;
end