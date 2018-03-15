%Based on AoA, AoD and ToA Characteristics of Scattered Multipath Clusters for THz Indoor Channel
%Modeling, Sebastian Priebe, Martin Jacob, Thomas Kurner 


num_channels=1;
num_clusters=1;

for k = 1:num_channels       % loop over number of channels
    % ToA and AoA of the reflectors (fixed values);
    ToA_ref= [10 20 30 40 50];%ns
    AoA_ref_azimuth=[10 20 30 40 50];%degree
    AoA_ref_elevation=[10 20 30 40 50];%degree
    
    
% parameters related to angle distribution of the scattering
a1_Az=0.583;
sigma1_Az=1.895;
x1_mu_Az=0;
a2_Az=0.417;
sigma2_Az=5.487;
x2_mu_Az=0;
data_points=100;
angle_start=-5;
angle_end=5;
angle_resolution=2;

a1_El=0.724;
sigma1_El=2.198;
x1_mu_El=0;
a2_El=0.276;
sigma2_El=7.297;
x2_mu_El=0;
data_points=100;
angle_start=-5;
angle_end=5;
angle_resolution=2;
% parameter relted to the roughness of the surface
A_mu=[0.246 0.346 0.123 0.564 0.800 ];% with roughness this value increase
A_max=1;
A_min=0.1;

% angle calculation
Azimuth=Angle(a1_Az,sigma1_Az,x1_mu_Az,a2_Az,sigma2_Az,x2_mu_Az,data_points,angle_start,angle_end,angle_resolution);
hist(Azimuth, 50);
Elavation=Angle(a1_El,sigma1_El,x1_mu_El,a2_El,sigma2_El,x2_mu_El,data_points,angle_start,angle_end,angle_resolution);
hist(Elavation, 50);

% relative reflection amplitude (represents roughness)
A_i=Relative_amp(A_mu,A_max,A_min)

%  reflection power

P_ref= 0.8;

% cluster power

Std_Azimuth = std(Azimuth);
Std_Elavation = std(Elavation);



% Time delay profile


T_specular=0;%in ps
% ToA is generated using GMM

a1_ToA_AZ=0.661;
sigma1_ToA_AZ=1.196;
x1_mu_ToA_AZ=.570;
a2_ToA_AZ=1.108;
sigma2_ToA_AZ=14.430;
x2_mu_ToA_AZ=-12.525;
data_points_AZ=200;
angle_start_AZ=0;
angle_end_AZ=30;
angle_resolution_AZ=1;

a1_ToA_EL=0.705;
sigma1_ToA_EL=1.018;
x1_mu_ToA_EL=0.928;
a2_ToA_EL=0.804;
sigma2_ToA_EL=11.594;
x2_mu_ToA_EL=-3.386;
data_points_EL=200;
angle_start_EL=0;
angle_end_EL=30;
angle_resolution_EL=1;

Azimuth_TOA=Angle(a1_ToA_AZ,sigma1_ToA_AZ,x1_mu_ToA_AZ,a2_ToA_AZ,sigma2_ToA_AZ,x2_mu_ToA_AZ,data_points_AZ,angle_start_AZ,angle_end_AZ,angle_resolution_AZ);
figure
hist(Azimuth_TOA, 50);
Elavation_TOA=Angle(a1_ToA_EL,sigma1_ToA_EL,x1_mu_ToA_EL,a2_ToA_EL,sigma2_ToA_EL,x2_mu_ToA_EL,data_points_EL,angle_start_EL,angle_end_EL,angle_resolution_EL);
figure 
hist(Elavation_TOA, 50);
A=[length(Azimuth_TOA) length(Elavation_TOA) length(Elavation) length(Azimuth)]
data_length=min(A)
TOA=zeros(1,data_length);


%Power of the cluster
P=zeros(1,data_length);
for i=1:length(Azimuth)
P(i)=P_ref*(A_i^2)*exp(-0.5*(Azimuth(i)^2/Std_Azimuth^2))*exp(-0.5*(Elavation(i)^2/Std_Elavation^2));
end

Elavation=zeros(1,data_length)
Azimuth_TOA=ones(1,data_length)
Elavation_TOA=ones(1,data_length)
for i=1:data_length
TOA(i)=T_specular+ Azimuth_TOA(i)*(Azimuth(i)^2)+(Elavation_TOA(i)*Elavation(i)^2);
end
figure 
stem(Azimuth,TOA)
tmp_t=TOA;
tmp_h=P.^0.5;

   np = data_length;  % number of rays (or paths) for this realization
  [sort_tmp_t,sort_ix] = sort(tmp_t(1:np));    % sort in ascending time order
  t(1:np) = sort_tmp_t;
  h(1:np) = tmp_h(sort_ix(1:np));
  
    figure 
    stem(t,h)
    
    
   ts=0.05; 
    
    h_ct=h;
t_ct=t;
num_channels=1
%[hN,N] = uwb_sv_cnvrt_ct_15_4a( h_ct, t_ct, np, num_channels, ts );
ts=ts;
min_Nfs = 100; % GHz
N = max( 1, ceil(min_Nfs*ts) );  % N*fs = N/ts is the intermediate sampling frequency before decimation
N = 2^nextpow2(N); % make N a power of 2 to facilitate efficient multi-stage decimation

Nfs = N / ts;
t_max = max(t(:));  % maximum time value across all channels
h_len = 1 + floor(t_max * Nfs);  % number of time samples at resolution ts / N
hN = zeros(h_len,num_channels);

  np_k = np;  % number of paths in this channel
  t_Nfs = 1 + floor(t(1:np_k) * Nfs);  % vector of quantized time indices for this channel
  for n = 1:np_k
    hN(t_Nfs(n)) = hN(t_Nfs(n)) + h_ct(n);
  end



if N > 1,
  h2 = resample(hN, 1, N);  % decimate the columns of hN by factor N
else
  h2 = hN;
end
t = [0:(length(h2)-1)] * ts; 
h2 = h2 * N;


kappa = 0;%1.53; 
  fc = 300;   % GHz
  fs = 20;   % 2 - 10 GHz
h=h2;
h_len = length(h(:,1))
f = [1-fs/fc/2 : fs/fc/h_len/2 : 1+fs/fc/2].^(-2*(kappa));
f = [f(h_len : 2*h_len), f(1 : h_len-1)]';
i = (-1)^(1/2);    % complex i
for  c = 1:num_channels
    % add the frequency dependency
   h3 = zeros(2*h_len, 1);
    h3(1 : h_len) = h(:,c);  % zero padding
    fh2 = fft(h3);
    fh2 = fh2 .* f;
   h3 = ifft(fh2);
    h(:,c) = h3(1:h_len);
    %}
    % change to complex baseband channel
    phi = rand(h_len, 1).*(2*pi);
    h(:,c) = h(:,c) .* exp(phi .* i);
    
    % Normalize the channel energy to 1
    h(:,c) = h(:,c)/sqrt(h(:,c)' * h(:,c) );
end


end