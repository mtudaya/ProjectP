% Saleh-Valenzuela model channel parameters
clear all;
channel_num=3;
[Lam_cluster,lambda_ray,Gam_cluster,gamma_ray,std_cluster_fading_db,std_ray_fading_db,nlos,std_shadow] = sv_channel_parameter(channel_num );

std_cluster = 1/sqrt(2*Lam_cluster); % std dev (nsec) of cluster arrival spacing
std_ray = 1/sqrt(2*lambda_ray); % std dev (nsec) of ray arrival spacing
std_cluster_fading = 10^(std_cluster_fading_db/20);%Standard deviation of log-normal variable for cluster fading (not db)
std_ray_fading = 10^(std_ray_fading_db/20);%Standard deviation of log-normal variable for ray fading

mu_const = (std_cluster_fading^2+std_ray_fading^2)*log(10)/20; % pre-compute
Nfs_min = 100e+9; % 100 GHz
fs= 6e+9;% sampling rate
ts=1/fs; % sampling time
ts_ns=ts*10^9; %sampling time in ns
NLOS = nlos;
Omega0=1;%mean energy of the first path of the first cluster
num_channels=1; %number of channels
%h_len =1000; % length of the channel
%h = zeros(h_len,num_channels);
%t = zeros(h_len,num_channels);
%t0 = zeros(1,num_channels);
%np = zeros(1,num_channels);

fprintf(1,['Model Parameters\n' ...
' Lam_cluster = %.4f, lambda_ray = %.4f, Gam_cluster = %.4f, gamma_ray = %.4f\n' ...
' std_cluster_fading_db = %.4f, std_ray_fading_db = %.4f, NLOS flag = %d, std_shadow = %.4f\n'], ...
Lam_cluster, lambda_ray, Gam_cluster, gamma_ray, std_cluster_fading_db, std_ray_fading_db, nlos, std_shadow);

for k = 1:num_channels % loop over number of channels
% calculation of arrival times
    if NLOS==0
        TC_cell = [0];% first cluster arrival time
        TC=0;
    else
        TC_cell = [randn^2/(2*Lam_cluster)+randn^2/(2*Lam_cluster)];% first cluster arrival time
        TC = (randn^2/(2*Lam_cluster)+randn^2/(2*Lam_cluster)); 
    end


while TC(end) <= 10*Gam_cluster
    dT = randn^2/(2*Lam_cluster)+randn^2/(2*Lam_cluster);
    TC = [TC TC(end)+dT];
end


for n = 1:length(TC)
   % tau(n) = {[T1(n)]};
    Tr(n) = {[TC(n)]}; % assign time of the cluster to time of first ray
    %Tr(n) = 0;
    while Tr{n}(end) <= Tr{n}(1)+10*gamma_ray %Tr{n}(end) <= Tr{n}(1)+10*gamma_ray
     
        dt = randn^2/(2*lambda_ray)+randn^2/(2*lambda_ray);
        Tr(n) = {[Tr{n} Tr{n}(end)+dt]};
    end
    %calculation of amplitude of reflection
    % for IEEE  P802.15-02 version
    
    mu(n) = {((10*log(Omega0))-(10*TC(n)/Gam_cluster)-(10*(Tr{n}-TC(n))/gamma_ray))/log(10) - mu_const};
    %
    beta(n) = {10.^((std_cluster_fading*randn(size(mu{n}))+std_ray_fading*randn(size(mu{n}))+mu{n})/20)};
    %alpha(n) = {beta{n}};  % 
    alpha(n) = {beta{n}.*sign(rand(size(beta{n}))-0.5)};  % Intel version, +/-1
    
end

    time_max = max(Tr{1});  % find maximum arrival time of rays
    for n = 2:length(TC)
        if max(Tr{n}) > time_max
            time_max = max(Tr{n});
        end
    end
% convert continuous time channel model TO discrete time

N = max( 1, ceil(Nfs_min*ts) ); % N*fs = N/ts is the intermediate sampling frequency before decimation
N = 2^nextpow2(N); % make N a power of 2 to facilitate efficient multi-stage decimation
Nfs = N / ts_ns;
h = zeros(1,floor(time_max*Nfs)+1);
hN = zeros(1,floor(time_max*Nfs)+1);
for n = 1:length(TC)
    TrN{n} = floor(Tr{n}*Nfs);  % quantized time indices
    %TrN{n} = 1+floor(Tr{n}*Nfs);  % quantized time indices
    
end
for n = 1:length(TC)
    for m = 1:length(Tr{n})
        %h(TrN{n}(m)+1-TrN{1}(1)) = h(TrN{n}(m)+1-TrN{1}(1))+alpha{n}(m);  % 
        hN(TrN{n}(m)+1-TrN{1}(1)) = hN(TrN{n}(m)+1-TrN{1}(1))+alpha{n}(m);  % 
    end
end

hN = N*resample(hN,1,N);
time_max = ceil((10*Gam_cluster+10*gamma_ray)/ts_ns);  % maximum arrival time for channel
    if length(hN) > time_max
        hN = hN(1:time_max);  % concatenate h to maximum channel length
    end
E = sum(hN.*conj(hN));  % compute total channel energy
hN = hN./sqrt(E);  % normalize total energy to 1

fac = 10^(std_shadow*randn/20);
hN = hN.*fac;
hN=abs(hN);
end
h_len = length(hN);
t = [0:(h_len-1)] * ts_ns; % compute time values
figure(1); clf; plot(t,hN); grid on
title('Impulse response realizations')
xlabel('Time (nS)')