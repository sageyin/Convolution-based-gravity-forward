function g = GraconvelP(Density,dr,r,t,Style)
%% Construct the circular kernel matrix and calculate the gravity field using FFT algorithm
% Editor：Xianzhe Yin 2022/9/05 China University of Geosciences(Beijing)
%% Parameters
% ===== input =====
% Density: Model density matrixr (unit:kg/m^3)
% dr: Relative distance between model and griddr (unit:m)
% r : Size of the sub-cell model in x,y,z direction, respectively (unit:m)
% t : The number of points in x,y,z direction,respectively 
% Style : Type of gravitational field, including the vertical component of gravity and its tensor
% ===== out =====
% g : the vertical component of gravity or its tensor

[sN,sW,sz] = size(Density);            % Source grid size   
dW= dr(1); dN = dr(2); dz = dr(3);     % Source grid size spacing（dx，dy，dz）unit:m 
rW = r(1);  rN = r(2);  rz = r(3);     % Relative distance of the observation grid to the first grid of the source（rx，ry，rz）unit:m  
tN = t(1);  tW = t(2);  tz = t(3);     % Observation grid size

%% ====== Size after mesh expansion = S + T - 1 ======
x = [0:tN-1,-sN+1:-1]*dN - rN;   
y = [0:tW-1,-sW+1:-1]*dW - rW;
z = [0:tz-1,-sz+1:-1]*dz - rz;            

%% ====== Construction of circular kernel matrix ======
[W,N,Z] = meshgrid(y,x,z);
K=Cal_tranGraf(N,W,Z,0,0,0,dN,dW,dz,1,Style);
F = fftn(K);                 

%% ====== Expand the source grid to the same scale as the circular kernel matrix ======
S = zeros(size(N));
S(1:sN,1:sW,1:sz) = Density;        % Source grid expansion 0

%% ====== Cyclic convolution operations are transformed into dot product operations in the frequency domain ======
T = ifftn(fftn(S).*F);
g = T(1:tN,1:tW,1:tz);      % Intercepting grid gravity anomaly data

end
