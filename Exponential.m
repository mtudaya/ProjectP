% generate some sample data for clustering model

clear all
num_channels=1;
L=10; % number of clusters
max_ray_percluster=30; % max ray per cluster
ts = 0.167; % sampling time (nsec)
%sigma_cluster = 1.75;% sigma_cluster  Standard deviation of normally distributed variable for cluster energy


alpha=zeros(L,1);
bete=zeros(L,1);
randn('state',12);    % initialize state of function for repeatability
rand('state',12);     % initialize state of function for repeatability


tmp_h = zeros(1000,1);
tmp_t = zeros(1000,1);
h_len=1000;
t_cluster=zeros(L,1);
h = zeros(h_len,num_channels);
t = zeros(h_len,num_channels);
t0 = zeros(1,num_channels);
np = zeros(1,num_channels);
cm_num = 6;  % channel model number from 1 to 8
% get channel model params based on this channel model number
[Lam,lambda,Lmean,lambda_mode,lambda_1,lambda_2,beta,Gam,gamma_0,Kgamma, ...
 sigma_cluster,nlos,gamma_rise,gamma_1,chi,m0,Km,sigma_m0,sigma_Km, ...
 sfading_mode,m0_sp,std_shdw,kappa,fc,fs] = uwb_sv_params_15_4a( cm_num );
fprintf(1,['Model Parameters\n' ...
  '  Lam = %.4f, lambda = %.4f, Lmean = %.4f, lambda_mode(FLAG) = %d\n' ...
  '  lambda_1 = %.4f, lambda_2 = %.4f, beta = %.4f\n' ...
  '  Gam = %.4f, gamma0 = %.4f, Kgamma = %.4f, sigma_cluster = %.4f\n' ...
  '  nlos(FLAG) = %d, gamma_rise = %.4f, gamma_1 = %.4f, chi = %.4f\n' ...
  '  m0 = %.4f, Km = %.4f, sigma_m0 = %.4f, sigma_Km = %.4f\n' ...
  '  sfading_mode(FLAG) = %d, m0_sp = %.4f, std_shdw = %.4f\n', ...
  '  kappa = %.4f, fc = %.4fGHz, fs = %.4fGHz\n'], ...
  Lam,lambda,Lmean,lambda_mode,lambda_1,lambda_2,beta,Gam,gamma_0,Kgamma, ...
  sigma_cluster,nlos,gamma_rise,gamma_1,chi,m0,Km,sigma_m0,sigma_Km,...
  sfading_mode,m0_sp,std_shdw,kappa,fc,fs);
std_lam = 1/sqrt(2*lambda); % std dev (nsec) of ray arrival spacing
std_L = 1/sqrt(2*Lam);      % std dev (nsec) of cluster arrival spacing
for k = 1:num_channels       % loop over number of channels
path_ix = 0;
cluster_ix = 0;
ray_ix = 0; 
Tc=0;
cluster_index = zeros(1,L);
for ncluster = 1:L
Tc = Tc + (std_L*randn)^2 + (std_L*randn)^2;
t_cluster(ncluster)=Tc;
gamma = Kgamma*Tc + gamma_0; % delay dependent cluster decay time
%gamma = gamma_0; % delay dependent cluster decay time
beta(ncluster)=gamma;
Mcluster = sigma_cluster*randn;
Pcluster = 10*log10(exp(-1*Tc/Gam))+Mcluster; % total cluster power
Pcluster = 10^(Pcluster*0.1);
alpha(ncluster)=Pcluster;
end

% Normalized cluster power
alpha = alpha / sqrt(alpha' * alpha);


for ncluster = 1:L
    Tr = 0;  % first ray arrival defined to be time 0 relative to cluster
    ray_ix = 0;
    cluster_index(ncluster) = path_ix+1; % remember the cluster location
    Tc=t_cluster(ncluster);
    cluster_ix = cluster_ix+1;
    while (Tr < 10*beta(ncluster) && ray_ix < max_ray_percluster),
      t_val = (Tc+Tr)  % time of arrival of this ray
      h_val = alpha(ncluster)*exp(-Tr/beta(ncluster)); %h_val = Pcluster/gamma*exp(-Tr/gamma);
      path_ix = path_ix + 1;  % row index of this ray
      ray_ix = ray_ix + 1;  % row index of this ray
      tmp_h(path_ix) = h_val;
      tmp_h_cluster(ray_ix,cluster_ix,k) = h_val;
      tmp_t(path_ix) = t_val;    
      tmp_t_cluster(ray_ix,cluster_ix,k) = t_val;    
      Tr = Tr + (std_lam*randn)^2 + (std_lam*randn)^2;
    end
   
    Tc = Tc + (std_L*randn)^2 + (std_L*randn)^2;
end
 np(k) = path_ix;  % number of rays (or paths) for this realization
 [sort_tmp_t,sort_ix] = sort(tmp_t(1:np(k)));    % sort in ascending time order
  t(1:np(k),k) = sort_tmp_t;
  h(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
end

h_ct=h;

ts = 1*1/fs;  % sampling frequency
% convert continuous-time channel model h_ct to N-times oversampled discrete-time samples
% h_ct, t, np, and num_channels are as specified in uwb_sv_model
% ts is the desired time resolution
% 
% hN will be produced with time resolution ts / N.
% It is up to the user to then apply any filtering and/or complex downconversion and then
% decimate by N to finally obtain an impulse response at time resolution ts.

min_Nfs = 100; % GHz
N = max( 1, ceil(min_Nfs*ts) );  % N*fs = N/ts is the intermediate sampling frequency before decimation
N = 2^nextpow2(N); % make N a power of 2 to facilitate efficient multi-stage decimation

Nfs = N / ts;
t_max = max(t(:));  % maximum time value across all channels
h_len = 1 + floor(t_max * Nfs);  % number of time samples at resolution ts / N
hN = zeros(h_len,num_channels);
for k = 1:num_channels
  np_k = np(k);  % number of paths in this channel
  t_Nfs = 1 + floor(t(1:np_k,k) * Nfs);  % vector of quantized time indices for this channel
  for n = 1:np_k
    hN(t_Nfs(n),k) = hN(t_Nfs(n),k) + h_ct(n,k);
  end
end

if N > 1,
  h = resample(hN, 1, N);  % decimate the columns of hN by factor N
else
  h = hN;
end

h_len=length(h);
t = [0:(h_len-1)] * ts; 
plot(t,10*log10(abs(h).^2));
%plot(t,h);