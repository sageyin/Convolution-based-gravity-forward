%% The main function of the calculation of the vertical component of the gravity field or its tensors in space domain
% Editor：Xianzhe Yin 2022/9/05 China University of Geosciences(Beijing)
clear
close all
clc
% We pick the NED (North-East-Downward) coordinate system as positive x, y, and z directions.
%% ====== Defining the observation grid ======
ObsGrid.dn=2;ObsGrid.de=2;ObsGrid.dz=2;
ObsGrid.Nmin=-100+ObsGrid.dn/2; ObsGrid.Nmax=100-ObsGrid.dn/2; % North-South
ObsGrid.Emin=-150+ObsGrid.de/2; ObsGrid.Emax=150-ObsGrid.de/2; % West-East
ObsGrid.zmin=0; ObsGrid.zmax=0;
ObsGrid.n=ObsGrid.Nmin:ObsGrid.dn:ObsGrid.Nmax;
ObsGrid.e=ObsGrid.Emin:ObsGrid.de:ObsGrid.Emax;
ObsGrid.z=ObsGrid.zmin:ObsGrid.dz:ObsGrid.zmax; 
[ObsGrid.E,ObsGrid.N,ObsGrid.Z]=meshgrid(ObsGrid.e,ObsGrid.n,ObsGrid.z);

%% ====== Rectangular model construction ======
SouceGrid.dn=ObsGrid.dn;SouceGrid.de=ObsGrid.de;SouceGrid.dz=ObsGrid.dz;
SouceGrid.Nmin=-100; SouceGrid.Nmax=100; % North-South
SouceGrid.Emin=-150; SouceGrid.Emax=150; % West-East
SouceGrid.zmin=0; SouceGrid.zmax=100;
SouceGrid.n=SouceGrid.Nmin+SouceGrid.dn/2:SouceGrid.dn:SouceGrid.Nmax-SouceGrid.dn/2;
SouceGrid.e=SouceGrid.Emin+SouceGrid.de/2:SouceGrid.de:SouceGrid.Emax-SouceGrid.de/2;
SouceGrid.z=SouceGrid.zmin+SouceGrid.dz/2:SouceGrid.dz:SouceGrid.zmax-SouceGrid.dz/2; 
[SouceGrid.E,SouceGrid.N,SouceGrid.Z]=meshgrid(SouceGrid.e,SouceGrid.n,SouceGrid.z);% Dissection of the underground grid
SouceGrid.density=zeros(size(SouceGrid.E));

logp=logical(SouceGrid.E<=80 & SouceGrid.E>=40 & SouceGrid.N<=80 & SouceGrid.N>=40 ...
                      & SouceGrid.Z<=20 & SouceGrid.Z>=10 );    % Rectangle 50*50*10 Center position: (0, 0, 35)
Souce.E=SouceGrid.E(logp);
Souce.N=SouceGrid.N(logp );
Souce.Z=SouceGrid.Z(logp);
Souce.density=1000; % unit:kg/m^3
SouceGrid.density(logp)=Souce.density;

%% ====== Forward modelling of gravity field by our method in space domain ======

tic
Style='gz';
Souce.Num=length(Souce.E);
gg=0;
for n=1:Souce.Num
    g=Cal_tranGraf(ObsGrid.N,ObsGrid.E,ObsGrid.Z,Souce.N(n),Souce.E(n),Souce.Z(n),SouceGrid.dn,SouceGrid.de,SouceGrid.dz,Souce.density,Style);
    gg=gg+g;
end
toc

%% ====== Unit Conversions ======
gg=gg*10^5;   %  m/s^2 converted to mGal
%gg=gg*10^6;  %  m/s^2 converted to g.u
%gg=gg*1e9;   %  s^(-2) to E

%% ====== Visualization ======
figure()
contourf(ObsGrid.e,ObsGrid.n,gg)
colormap('jet');colorbar
axis equal
xlabel('West-East(m)');
ylabel('South-North(m)')
title('Gravity field calculated in space domain')

%% ====== Data Storage ======
save('mode01',"gg"); 
