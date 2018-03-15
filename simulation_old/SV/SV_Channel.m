% Saleh-Valenzuela model channel parameters



Lam_cluster = 0.2; %Cluster arrival rate (clusters per nsec)
std_cluster = 1/sqrt(2*Lam_cluster); % std dev (nsec) of cluster arrival spacing
lambda_ray = 0.3;% Ray arrival rate (rays per nsec)
std_ray = 1/sqrt(2*lambda_ray); % std dev (nsec) of ray arrival spacing
Gam_cluster = 5.5;%Cluster decay factor (time constant, nsec)
gamma_ray = 6.7;% Ray decay factor (time constant, nsec)
std_cluster_fading_db =3.3941;%Standard deviation of log-normal variable for cluster fading 
std_ray_fading_db=3.3941;%Standard deviation of log-normal variable for ray fading
std_cluster_fading = 10^(std_cluster_fading_db/20);%Standard deviation of log-normal variable for cluster fading (not db)
std_ray_fading = 10^(std_ray_fading_db/20);%Standard deviation of log-normal variable for ray fading
std_shadow = 3;%Standard deviation of log-normal shadowing of entire impulse response
mu_const = (std_cluster_fading^2+std_ray_fading^2)*log(10)/20; % pre-compute
Nfs_min = 100e+9; % 100 GHz
fs= 6e+9;% sampling rate
ts=1/fs; % sampling time
ts_ns=ts*10^9; %sampling time in ns
LOS = 1;
Omega0=1;%mean energy of the first path of the first cluster
num_channels=1; %number of channels
%h_len =1000; % length of the channel
%h = zeros(h_len,num_channels);
%t = zeros(h_len,num_channels);
%t0 = zeros(1,num_channels);
%np = zeros(1,num_channels);

for k = 1:num_channels % loop over number of channels
% calculation of arrival times
    if LOS==0
        TC = [0];% first cluster arrival time
    else
        TC = [randn^2/(2*Lam_cluster)+randn^2/(2*Lam_cluster)];% first cluster arrival time
        
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
hN = hN(1:time_max);  % concatenate h to maximum channel length

E = sum(hN.*conj(hN));  % compute total channel energy
hN = hN./sqrt(E);  % normalize total energy to 1

fac = 10^(std_shadow*randn/20);
hN = hN.*fac;
hN=abs(hN);
end
E=hN*(hN)';
h_len = length(hN);
t = [0:(h_len-1)] * ts_ns; % compute time values
figure(1); clf; plot(t,hN); grid on
title('Impulse response realizations')
xlabel('Time (nS)')