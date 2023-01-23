%% The main function of the new algorithm, the calculation of the vertical component of the gravity field or its tensors
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
                      & SouceGrid.Z<=20 & SouceGrid.Z>=10 );    % Rectangle 50*50*20 Center position: (0, 0, 35)
Souce.E=SouceGrid.E(logp);
Souce.N=SouceGrid.N(logp );
Souce.Z=SouceGrid.Z(logp);
Souce.density=1000; % unit:kg/m^3
SouceGrid.density(logp)=Souce.density;

%% ====== Forward modelling of gravity field by our method ======
tic
dr=[SouceGrid.de,SouceGrid.dn,SouceGrid.dz];
r=[0,0,SouceGrid.dz/2]; % Upward is positive
t=[0,0,1];
t(1)=size(ObsGrid.E,1);t(2)=size(ObsGrid.E,2);
g =GraconvelP(SouceGrid.density,dr,r,t,'gz');
toc

%% ====== Unit Conversions ======
g=g*10^5;  %   m/s^2 converted to mGal
%g=g*10^6; %   m/s^2 converted to g.u
%g=g*1e9;  %   s^(-2) to E

%% ====== load data from spatial domain calculations ======
g0=load('mode01',"gg");  % Data obtained from spatial domain calculations
er=g0.gg-g;              % Calculation of the difference between the two methods

%% ====== Visualization ======
figure()
contourf(ObsGrid.e,ObsGrid.n,g)
colormap('jet');colorbar
axis equal
xlabel('West-East(m)');
ylabel('South-North(m)')
title('Gravity field calculated by our method')

figure()
contourf(ObsGrid.e,ObsGrid.n,er,50,'LineStyle','none');
colormap('jet');colorbar
axis equal
xlabel('West-East(m)');
ylabel('South-North(m)')
title('The difference between the two methods')

figure()
histogram(er)
xlabel('misfit(mGal)');
ylabel('points');
title('Deviation statistics')




