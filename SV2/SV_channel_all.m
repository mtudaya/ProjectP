
%[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=SV_channel_all(num_clusters_fc);
function [h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=SV_channel_all(num_clusters_fc);

Lam = 0.8;  lambda = 10;  Lmean = 2;
  lambda_mode = 2;  %0 -> Poisson process for the ray arrival times, 2 channel tap
  lambda_1 =  NaN; lambda_2 =  NaN; beta = NaN; 
  % MPC decay
  Gam =10; % higher gam low attenuation
  gamma_0 =  1.5; Kgamma = 0; sigma_cluster = 2.93;
  nlos = 1; 
  gamma_rise = NaN; gamma_1 = NaN; chi = NaN; % dummy in this scenario
  % Small-scale Fading
 % m0 = 0.69; Km = 0; sigma_m0 = 0.32; sigma_Km = 0;
 % sfading_mode = 0; m0_sp = NaN;
  m0 = 10; Km = 0; sigma_m0 = 10; sigma_Km = 0;
  sfading_mode = 0; m0_sp = NaN;
  % Large-scale Fading -- Shadowing
  std_shdw = 3 %3.51;
  % Frequency Dependence
  kappa = 0;%1.53; 
  fc = 6;   % GHz
  fs = 8;   % 2 - 10 GHz
  ts=1/fs;
num_channels=1;
%num_clusters_fc=4;
% initialize and precompute some things
std_L = 1/sqrt(2*Lam);      % std dev (nsec) of cluster arrival spacing
std_lam = 1/sqrt(2*lambda); % std dev (nsec) of ray arrival spacing
h_len = 1000;  % there must be a better estimate of # of paths than this
ngrow = 1000; % amount to grow data structure if more paths are needed
h = zeros(h_len,num_channels);
t = zeros(h_len,num_channels);
t0 = zeros(1,num_channels);
np = zeros(1,num_channels);
cluster_non_overlap=0; % if clusteres are not overlap 1

for k = 1:num_channels       % loop over number of channels
   
  tmp_h = zeros(size(h,1),1);
  tmp_t = zeros(size(h,1),1);
  if nlos == 1,
    Tc = (std_L*randn)^2 + (std_L*randn)^2;  +1;  % First cluster random arrival
    % make Tc as multiple of sampling time
    Tc= (ceil(Tc/ts)-1)*ts;
  else
    Tc = 0;                   % First cluster arrival occurs at time 0
  end
  t0(k) = Tc;
  
  
  L=num_clusters_fc;
  disp(L);  
  cluster_index = zeros(1,L);
  path_ix = 0;
  
  nak_m = [];

  
  for ncluster = 1:L
      
    % Determine Ray arrivals for each cluster
    Tr = 0;  % first ray arrival defined to be time 0 relative to cluster
    actual_cluster_times(ncluster,k)=Tc;
    cluster_index(ncluster) = path_ix+1 ; % remember the cluster location
    ray_ix = 0;
    gamma = Kgamma*Tc + gamma_0; % delay dependent cluster decay time
    if nlos == 2 & ncluster == 1
      gamma = gamma_1;
    end
    Mcluster = sigma_cluster;%*randn;
    Pcluster = 10*log10(exp(-1*Tc/Gam));%+Mcluster; % total cluster power
    Pcluster = 10^(Pcluster*0.1)+0.8;
    Pcluster_actual(ncluster,k)=Pcluster; 
    if ncluster ~=1 
        Pcluster = Pcluster_actual(ncluster-1,k)/15;
        Pcluster_actual(ncluster,k)=Pcluster; 
    end
    while (Tr < 15*gamma),
      t_val = (Tc+Tr);  % time of arrival of this ray
      ray_ix=ray_ix+1;
      actual_ray_times(ray_ix,ncluster,k)=t_val;
     %h_val = Pcluster/gamma*exp(-Tr/gamma);
      h_val = Pcluster*exp(-Tr/gamma)
      actual_ray_times_power(ray_ix,ncluster,k)=h_val;
     % end
      path_ix = path_ix + 1;  % row index of this ray
      if path_ix > h_len,
        % grow the output structures to handle more paths as needed
        tmp_h = [tmp_h; zeros(ngrow,1)];
        tmp_t = [tmp_t; zeros(ngrow,1)];
        h = [h; zeros(ngrow,num_channels)];
        t = [t; zeros(ngrow,num_channels)];
        h_len = h_len + ngrow;
      end
      tmp_h(path_ix) = h_val;
      tmp_t(path_ix) = t_val;
      if lambda_mode == 0
        Tr = Tr + (std_lam*randn)^2 + (std_lam*randn)^2;
        %Tr= (ceil(Tr/ts)-1)*ts;
      elseif lambda_mode == 2
        Tr = Tr + ts;
      else
        error('lambda mode is wrong!')
      end
      % generate log-normal distributed nakagami m-factor
      m_mu = m0 - Km*t_val;
      m_std = sigma_m0 - sigma_Km*t_val;
      nak_m = [nak_m, lognrnd(m_mu, m_std)];
    end
    
    
    if cluster_non_overlap==1
        while 1
            Tc = Tc + (std_L*randn)^2 + (std_L*randn)^2;
            if Tc > t_val
            break
            end
        end
    else
        %Tc = Tc + (std_L*randn)^2 + (std_L*randn)^2;
        while 1
            Tc_old=Tc;
            Tc = Tc + 3 %(std_L*randn)^2 + (std_L*randn)^2;
            if Tc > Tc_old+2;
                break
            end
        end
    end
    
  end

  % change m value of the first multipath to be the deterministic value
  if sfading_mode == 1
    nak_ms(cluster_index(1)) = m0_sp;
  elseif sfading_mode == 2
    nak_ms(cluster_index) = m0_sp;
  end
  
%   apply nakagami
  for path = 1:path_ix
     h_val = (gamrnd(nak_m(path), tmp_h(path)/nak_m(path))).^(1/2);
     tmp_h(path) = h_val;
  end

  np(k) = path_ix;  % number of rays (or paths) for this realization
  [sort_tmp_t,sort_ix] = sort(tmp_t(1:np(k)));    % sort in ascending time order
  t(1:np(k),k) = sort_tmp_t;
  h(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
  h_normalized(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
  % now impose a log-normal shadowing on this realization
  %fac = 10^(std_shdw*randn/20) / sqrt( h(1:np(k),k)' * h(1:np(k),k) );
  %h(1:np(k),k) = h(1:np(k),k) * fac;
 % Normalize the channel energy to 1
   h_normalized(:,k) = h_normalized(:,k)/sqrt(h_normalized(:,k)' * h_normalized(:,k) );
  % h(:,k) = h(:,k)/sqrt(h(:,k)' * h(:,k) );
   %h_normalized2=h_normalized(:,k);
    %actual_cluster_times2=actual_cluster_times(:,k);
    %[sharedVals,idxsIntoA] = intersect(t(:,k),actual_cluster_times2);
     %Pcluster_actual(:,k)=h_normalized2(idxsIntoA); 
    %Pcluster_actual(:,k)= Pcluster_actual(:,k).^2;
   
end
%h(:,k) = h(:,k)/sqrt(h(:,k)' * h(:,k));
h_ct=h;
t_ct=t;

%[hN,N] = uwb_sv_cnvrt_ct_15_4a( h_ct, t_ct, np, num_channels, ts );
ts=ts;
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
  h2 = resample(hN, 1, N);  % decimate the columns of hN by factor N
else
  h2 = hN;
end
t = [0:(length(h2)-1)] * ts; 
h2 = h2 * N;



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
figure(1);
subplot 221
stem(t_ct,h_ct,'filled','markersize',3)
grid on
xlabel 'Time',ylabel 'Before Sampling'
subplot 222
stem(t,h2,'filled','markersize',3)
grid on
xlabel 'Time',ylabel 'Decimated'
subplot 223
stem(t,abs(h).^2,'filled','markersize',3)
grid on
xlabel 'Time',ylabel 'baseband'

subplot 224
plot(t,10*log10(abs(h).^2))
axis([0 inf -100 10])
grid on
xlabel 'Time',ylabel 'baseband'

figure(2);
subplot 211
stem(actual_ray_times,actual_ray_times_power,'filled','markersize',3)
grid on
xlabel 'Time',ylabel 'Channel Impulse Response'
subplot 212

stem(actual_ray_times,10*log10(actual_ray_times_power),'filled','markersize',3)
grid on
xlabel 'Time',ylabel 'Channel Power DB'
figure(3);
plot(actual_ray_times,10*log10(actual_ray_times_power))
grid on
xlabel 'Time',ylabel 'Channel Power DB'

return