%% 计算laplace反演
%function [PSIO] = laplace_transform(forcing,SC,CRI,dx,dy,p,lat,f,type)

% 输入：
%forcing 强迫项
%SC 超松弛迭代因子，一般取1~2
%CRI 迭代判定条件，一般取10^-12
%dx,dy 水平方向经度（纬度）格点差
%p 输入垂直方向插值前气压场
%lat 纬度
%f 科式参数
%type =1 输入未处理的热力强迫， type=2, 输入动力强迫

function [PSI0,RR] = laplace_transform(forcing,SC,CRI,dx,dy,p,lat,f,type)
[nx,ny,np]=size(forcing);%dims for Input forcing
PSI0=zeros(nx,ny,np);
pp=(100000.0:-5000.0:10000.0);%插值后垂直方向气压层19层
dp=-5000.0;
nz=length(pp);
clear t1 t2;
%new_f=repmat(f,[nx,1,nz]);
new_f(1:nx,1:ny,1:nz)=2*7.27/10^5*sqrt(2)/2;
new_lat=repmat(lat,[nx,1,nz]);

rate=[1.497942e-06,1.727951e-06,1.799234e-06,2.371034e-06,2.346085e-06,2.781826e-06,4.076014e-06,9.296354e-06,1.807938e-05,3.234908e-05,5.818076e-05,0.0001368937];
R1=spline(p,rate,pp);
R=1./R1;
Ra=6.4E+06;

new_R1=zeros(nx,ny,nz);
new_pp=zeros(nx,ny,nz);

for iz=1:nz
    new_R1(:,:,iz)=R1(1,iz);
    new_pp(:,:,iz)=pp(1,iz);
end;

FRC0=zeros(nx,ny,nz);
for ix=1:nx
    for iy=1:ny
        FRC0(ix,iy,:)=spline(p,forcing(ix,iy,:),pp);
    end;
end;
           % t_f=nanmean(FRC0,1);
           % figure(2);
           % contourf(squeeze(t_f)');
           
if (type==1)
    s1=FRC0./(new_pp.*new_R1);
    [f_dx,f_dy,f_dp]=gradient(s1);
    s2=f_dp/dp;
    FRC1=-new_f.*s2*287;
    clear s1 s2;
end;
if (type==2)
    FRC1=FRC0;
end;
FRC=FRC1.*new_f*Ra*Ra.*cos(new_lat).*cos(new_lat);
            %t_f=nanmean(FRC,1);
            %figure(3);
            %contourf(squeeze(t_f)');
%PS0=FRC0*1E+5;;
PS0=zeros(nx,ny,nz);
PS1=zeros(nx,ny,nz);
for k1=1:100
    PS1=PS0;
    for iz=2:nz-1
        for iy=2:ny-1
            for ix=2:nx-1
               %%{
            COE011=1./(dx*dx);
            COE211=1./(dx*dx);
            COE101=cos(new_lat(ix,iy,iz))*sin(new_lat(ix,iy,iz))/(2*dy)+cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))/(dy*dy);
            COE121=-cos(new_lat(ix,iy,iz))*sin(new_lat(ix,iy,iz))/(2*dy)+cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))/(dy*dy);
            COE110=-Ra*Ra*(cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(ix,iy,iz))*(R(iz+1)-R(iz-1))/(2*dp)/(2*dp)+Ra*Ra*cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(iz,iy,iz)*R(iz)/(dp*dp);
            COE112=Ra*Ra*(cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(ix,iy,iz))*(R(iz+1)-R(iz-1))/(2*dp)/(2*dp)+Ra*Ra*cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(iz,iy,iz)*R(iz)/(dp*dp);
            COE111=-2/(dx*dx)-2*cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))/(dx*dx)-2*Ra*Ra*cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(ix,iy,iz)*R(iz)/(dp*dp);
            RR=FRC(ix,iy,iz)-(COE011*PS0(ix-1,iy,iz)+COE111*PS0(ix,iy,iz)+COE211*PS0(ix+1,iy,iz)+COE101*PS0(ix,iy-1,iz)+COE121*PS0(ix,iy+1,iz)+COE110*PS0(ix,iy,iz-1)+COE112*PS0(ix,iy,iz+1));
                %{       
                COE011=1./(dx*dx);
                COE111=-0.5*(4*COE011+(cos(new_lat(ix,iy,iz))*COE011)*(cos(new_lat(ix,iy-1,iz))+2*cos(new_lat(ix,iy,iz))+cos(new_lat(ix,iy+1,iz)))+Ra*Ra*cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(ix,iy,iz))/(dp*dp)*(R(iz-1)+2*R(iz)+R(iz+1));
                COE211=COE011;
                COE101=0.5*(cos(new_lat(ix,iy,iz))*COE011)*(cos(new_lat(ix,iy-1,iz))+cos(new_lat(ix,iy,iz)));
                COE121=0.5*(cos(new_lat(ix,iy,iz))*COE011)*(cos(new_lat(ix,iy,iz))+cos(new_lat(ix,iy+1,iz)));
                COE110=0.5*(Ra*Ra*cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(ix,iy,iz))/(dp*dp)*(R(iz-1)+R(iz));
                COE112=0.5*(Ra*Ra*cos(new_lat(ix,iy,iz))*cos(new_lat(ix,iy,iz))*new_f(ix,iy,iz)*new_f(ix,iy,iz))/(dp*dp)*(R(iz)+R(iz+1));
                RR=FRC(ix,iy,iz)-(COE011*PS0(ix-1,iy,iz)+COE111*PS0(ix,iy,iz)+COE211*PS0(ix+1,iy,iz)+COE101*PS0(ix,iy-1,iz)+COE121*PS0(ix,iy+1,iz)+COE110*PS0(ix,iy,iz-1)+COE112*PS0(ix,iy,iz+1));
               %}
 PS0(ix,iy,iz)=PS0(ix,iy,iz)+SC/COE111*RR;
            end;
        end;
    end;

    if (type==2)
        PS0(:,:,1)=PS0(:,:,2);
        PS0(:,:,nz)=PS0(:,:,nz-1);
    end;

    if (type==1)
        PS0(:,:,1)=PS0(:,:,2)+287*dp/pp(1)*FRC0(:,:,1);
        PS0(:,:,nz)=PS0(:,:,nz-1)-287*dp/pp(1)*FRC0(:,:,nz);
    end;
    if (all((abs(PS0-PS1))<=CRI))
        break;
    end;
end;
k1
PSI=PS1;

for iy=1:ny
    for ix=1:nx
      PSI0(ix,iy,:)=spline(pp,PSI(ix,iy,:),p);
    end;
end;


            
