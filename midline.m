clear all;
cd WT;
foldername='005';
load (sprintf('%s/velocityField_2D_%s.mat',foldername,foldername));
load(sprintf('%s/g1_%s.mat',foldername,foldername));
load(sprintf('%s/analysis.mat',foldername))
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
d_p_d_z_set(:,1:end-1)=(fp(:,2:end)-fp(:,1:end-1))/dz;time=size(u_original,1);
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
filename='xy_cy1p_005_tp019.csv';
t_start=19
midline_coord=readmatrix(sprintf('%s/%s',foldername,filename));

%%
for (t=t_start:size(u_original,1))
    im=imread(sprintf('%s/cylinder1_proper_%s/cmp_1_1_T%04d.tif',foldername,foldername,t_start));
    imshow (im,[0 max(max(im))]);
    hold on;
    plot (q,p,'linewidth',2,'color','r');
for (j=1:size(p,1))
    z_act(j,t)=interp1(zp(1,:),z_range,q(j)-bdry{1}(1),'spline',0/0);
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
    f_z(:,1)=fp(:,floor (z_act(j)-z_range(1)+1))*(1-(z_act(j)-floor (z_act(j))))+fp(:,floor (z_act(j))-z_range(1)+2)*(z_act(j)-floor (z_act(j)));
    gff_int(:,1)=gff(:,floor (z_act(j)-z_range(1)+1))*(1-(z_act(j)-floor (z_act(j))))+fp(:,floor (z_act(j))-z_range(1)+2)*(z_act(j)-floor (z_act(j)));
    phi(j,t)=interp1(f_z,phi_range,p(j),'spline',0/0);
end
   v_p=v_original{t};
   v_q=u_original{t};
  
   
end
%%
