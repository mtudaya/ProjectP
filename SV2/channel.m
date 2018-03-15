clear all;
randn('state',12);    % initialize state of function for repeatability
rand('state',12);     % initialize state of function for repeatability
num_channels=1;
nlos=1;
cluster_non_overlap=0; % if clusteres are not overlap value= 1 else 0
Lam= 0.8;%   Cluster arrival rate (clusters per nsec)
lambda=1;% Ray arrival rate (rays per nsec)
Lcls=2;% number of Clusters
Gam=30;%    Cluster decay factor (time constant, nsec)
gamma=4; %Ray decay factor (time constant, nsec)
m0 = 0.69; Km = 0; sigma_m0 = 0.32; sigma_Km = 0;
sfading_mode = 0; m0_sp = NaN
sigma_cluster = 2.93;%Standard deviation of normally distributed variable for cluster energy
lambda_mode = 1 % Poison arrival process, if it is =1 then tap delayed model
fs = 8;   % 2 - 10 GHz
ts = 1/fs;  % sampling frequency
 nak_m = []; 
for k = 1:num_channels  

std_L = 1/sqrt(2*Lam);      % std dev (nsec) of cluster arrival spacing
std_lam = 1/sqrt(2*lambda); % std dev (nsec) of ray arrival spacing
h_len = 4000; 

h = zeros(h_len,1);
t = zeros(h_len,1);
tmp_h = zeros(size(h,1),1);
tmp_h2=zeros(size(h,1),1);
tmp_t = zeros(size(h,1),1);
Tc = 0;  % First cluster random arrival
t0(k) = Tc;
L=Lcls;
disp(L);  
cluster_index = zeros(1,L);
path_ix = 0;
time_id=0;
for ncluster = 1:L
      
    % Determine Ray arrivals for each cluster
    Tr = 0;  % first ray arrival defined to be time 0 relative to cluster
    actual_cluster_times(ncluster,k)=Tc;
    cluster_index(ncluster) = path_ix+1 ; % remember the cluster location
    ray_ix = 0;
    Mcluster = sigma_cluster*randn;
    Pcluster = 10*log10(exp(-1*Tc/Gam));%+Mcluster; % total cluster power
    Pcluster = 10^(Pcluster*0.1)
    Pcluster_actual(ncluster,k)=Pcluster; 
    
    while (Tr < 5*gamma),
      t_val = (Tc+Tr);  % time of arrival of this ray
      ray_ix=ray_ix+1;
      actual_ray_times(ray_ix,ncluster,k)=t_val;
      h_val = Pcluster*exp(-Tr/gamma);
      actual_ray_h(ray_ix,ncluster,k)=h_val;
      path_ix = path_ix + 1  % row index of this ray
      tmp_h(path_ix) = h_val;
      time_id=t_val/ts;
      time_id=time_id+1;
      tmp_h2(time_id) = tmp_h2(time_id)+h_val;
      tmp_t(path_ix) = t_val;
     % time_id=time_id+1;
      if lambda_mode == 0
        Tr = Tr + (std_lam*randn)^2 + (std_lam*randn)^2;
      
      elseif lambda_mode == 1
        Tr = Tr + ts
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
        while 1
            Tc_old=Tc;
            Tc = Tc + ts*50;%(std_L*randn)^2 + (std_L*randn)^2;
            Tc=ceil(Tc/ts)*ts;
            if Tc > Tc_old+4;
            break
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
  for path = 1:time_id
     h_val = (gamrnd(nak_m(path), tmp_h2(path)/nak_m(path))).^(1/2);
    tmp_h2(path) = h_val;
  end

  np(k) = time_id;  % number of rays (or paths) for this realization
  [sort_tmp_t,sort_ix] = sort(tmp_t(1:np(k)));    % sort in ascending time order
  t(1:np(k),k) = sort_tmp_t;
  h(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
% Normalize the channel energy to 1
  
   %plot(actual_ray_times(:,ncluster),actual_ray_h(:,ncluster));hold on
    

   
end
end


h2=h;


h2(h2==0) = []; % remove zeros
time=(0:1:length(h2)-1)*ts;
t=time;

h_len=length(h2);
 % change to complex baseband channel
i = (-1)^(1/2);    % complex i
phi = rand(h_len, 1).*(2*pi);
h2 = h2 .* exp(phi .* i);
 
plot(time,10*log(abs(h2.^2)));
 
%{   

r=h2;

conv_h=zeros(2*(length(r))-1,length(r));

for i=1:length(r)
    conv_h(1:length(r),i)=r;
    conv_h(:,i)=circshift(conv_h(:,i),i-1);
end
H=conv_h;
Hm=conv_h';
H_Hm=Hm*H;

e= eig(H_Hm);
h=h2;
t=time;
for no_of_channels=1:num_channels
ray_index_array_cls_all=[];% estimated index of the cluster arrival  with different minimum time gap of the cluster arrivals
p_est_all_cls =[];% estimated power for all clusters 
h_data_all=(abs(h).^2);% calculation of actual channel power for all channels
h_data=(abs(h(:,no_of_channels)).^2);% calculation of actual channel power
min_time=3;%expected minimum time between clusters
%[h_data2]=h_data_smooth(h_data);
%h_data=h_data2;
h_data(h_data==0) = NaN; 
h_data_log=log(h_data);
min_h_data_log=min(h_data_log);
h_data_log(isnan(h_data)) = min_h_data_log; % adjusting zero power values
h_len = size(h_data,1);
% genrate measuremnent matrix
n_l=length(h_data_log); % length of the data
%L_data= length (h_data_log);%length of the channel
tr_expected=ts;% time resolution of the rays ( assume it is same as sampling rate)
min_sample= round((min_time/tr_expected),0); % minimum number of sample gap between two clusters
% select the data point which are local maximum
% assuming there is one cluster and claculate the alpha and beta (one cluster model).
h_data_mod_log=h_data_log;


Y=h_data_log;
const=1.0;
[alpha_all,beta_all,time_w,w_L2]=L2_opt(Y,tr_expected,0);
% generate the data using alpha and beta
p_est_all=zeros(n_l,1);
for i=1:length(time_w)
   p_est_all(i)= alpha_all*exp(const*beta_all*time_w(i));
end

%finding the cluster starting point (if the actual power is higher than the
%power of the one cluster model assume it is a cluster starting point

ray_index_array_final=[];
%finding the cluster starting point (if the actual power is higher than the
%power of the one cluster model assume it is a cluster starting point
ray_index=1;
num_iter=1;
while ray_index+min_sample<length(h_data)

if num_iter==1
    ray_index=1;
  
end


j=1;
ray_index_array=[];
tmp_h_pwd=[];
while ray_index<length(h_data)
    
while true && ray_index<length(h_data)
   if h_data(ray_index)> p_est_all(ray_index);
       
   break ; 
   else
   ray_index=ray_index+1;
   end
end

ray_index_array(j)=ray_index% cluster starting point
tmp_h_pwd(j)=h_data(ray_index);
%{
ray_index_array_cls_all(j,m)=ray_index;
ray_index_array_cls_channel_all(j,m,no_of_channels)=ray_index;
cluster_arrival_times(j,m,no_of_channels)= ray_index* (tr_expected);
%}
j=j+1;
ray_index=ray_index+min_sample; % shift ray index as a length of minimum cluster time gap
end

if length (ray_index_array)==1
    ray_index_array_final(num_iter+1)=ray_index_array(1);
    break
end

if num_iter==1
    ray_index_array_final(1)=ray_index_array(1);
end
% find the channel amplitude at the cluster starting points
%tmp_h_pwd=zeros(length(ray_index_array),1);
%tmp_h_pwd=h_data(ray_index_array);

%ray_index array contains the cluster starting points, based on these
%points we are finding next cluster.
y_data_temp=[];
alpha_temp=[];
beta_temp=[];
error_temp=[]
for arrayid=2:length(ray_index_array)
    % calculate the fitting error ofeach cluster compared to first cluster
    y_data_temp=h_data_log(ray_index_array(1):ray_index_array(arrayid)-1) ;
    [alpha_temp(arrayid-1),beta_temp(arrayid-1),time_w,weight_L2]=L2_opt(y_data_temp,tr_expected,0);  
    power_temp=alpha_temp(arrayid-1)*exp((beta_temp(arrayid-1).*time_w));
    error_temp(arrayid-1)=sum(weight_L2.*(y_data_temp-power_temp).^2)/length(weight_L2);
end

 [sort_error_temp, sort_error_temp_id] = sort(error_temp)   % sort 
 disp(ray_index_array);
 ray_index_array_final(num_iter+1)=ray_index_array(sort_error_temp_id(1)+1)
 num_iter=1+num_iter;
 ray_index= ray_index_array_final(num_iter)+min_sample
end


ray_index_array=[];
ray_index_array=ray_index_array_final;







for m=1:length(ray_index_array)
%calculate alpha and beta for clusters
alpha=zeros(length(m),1);
beta=zeros(length(m),1);

p_est_all_cluster=zeros(length(h_data),m);
Qc_indi=zeros(length(h_data),m);
error_cls=[];
%p_est_all_cls=zeros(length(h_data),1);
%time_w_ovp = zeros(length(h_data),length(ray_index_array));

%calculate alpha and beta for clusters
for cluster_index=1:m
     
  if  cluster_index==m % check this is the last cluster, if last cluster get data from start of the cluster to end of h_data
        Y=h_data_log(ray_index_array(cluster_index):length(h_data)) ;    
  else
  
        Y=h_data_log(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1);
  end
  [alpha(cluster_index),beta(cluster_index),time_w,w_L2]=L2_opt(Y,tr_expected,0);  
  alpha_cls_all(cluster_index,m,no_of_channels)=alpha(cluster_index);
  beta_cls_all(cluster_index,m,no_of_channels)=beta(cluster_index); 
  w(cluster_index,m)=length(w_L2);
 % calculate the power for a clusters   
    for k=1:length(time_w)
        p_est_all_cluster((k+ray_index_array(cluster_index)-1),cluster_index)= alpha(cluster_index)*exp(beta(cluster_index)*time_w(k));
        Qc_indi((k+ray_index_array(cluster_index)-1),cluster_index)= h_data(k);
        %calculate average error for this cluster#
        
       % error_cls((k+ray_index_array(cluster_index)-1),cluster_index)=w_L2(k)*(((h_data(k+ray_index_array(cluster_index)-1))-(alpha(cluster_index)*exp(beta(cluster_index)*time_w(k))))^2);
        error_cls(k,cluster_index)=w_L2(k)*(((h_data(k+ray_index_array(cluster_index)-1))-(alpha(cluster_index)*exp(beta(cluster_index)*time_w(k))))^2);
    end
    %calculate average error for this cluster#
    error_cls_sum(cluster_index,m)=sum(error_cls(:,cluster_index));
    Gc(cluster_index)=(sum(p_est_all_cluster(:,cluster_index)));
    Qc(cluster_index)=(sum(Qc_indi(:,cluster_index)));
    
end

p_est_all_cls(:,m) = sum(p_est_all_cluster,2); % getting the estimated power for all estimated clusters
p_est_all_cls_channel(:,m,no_of_channels) = sum(p_est_all_cluster,2);
p_est_all_channel(:,no_of_channels)=p_est_all; % estimated power of assuming one cluster
clearvars p_est_all_cluster;
Gc=log(Gc);
Qc=log(Qc);
Qc_mean=mean(Qc);
R(m)=1-((sum((Qc-Gc).^2))/(sum((Qc-Qc_mean).^2)))


%p_est_all_log=log(p_est_all); 
%p_est_all_cls_log=log(p_est_all_cls(:,m));
%no_of_clusters(m)=nnz(ray_index_array_cls_all(:,m)); % calculate number of clusters for all time gap values
%num_of_est_clusters(m,no_of_channels)=(nnz(ray_index_array_cls_all(:,m)))';% calculate number of clusters for all time gap values for all channels
h_data(isnan(h_data)) = 0;

min_h_data=10^-15; % filter out data less than -150 dB
RMSE_weight=ones(length(h_data),1);
for weightindex=1:length(h_data)
    if h_data(weightindex)>min_h_data;
       RMSE_weight(weightindex)=1;
    else
       RMSE_weight(weightindex)=0; 
    end
end



error=(h_data-p_est_all_cls(:,m));  % Errors
errorwwww(:,m,no_of_channels)=error;
error_squ=error.^2 ;  % Squared Error
error_squ = nonzeros(error_squ);% remove zero values of the error
errorwwww2(:,m,no_of_channels)=error_squ;
error_mean=mean(error_squ) ;  % Mean Squared Error
RMSE(m,no_of_channels) = sqrt(error_mean)  % Root Mean Squared Error

error2=RMSE_weight.*(h_data-p_est_all_cls(:,m)); 
error_squ2=error2.^2 ;  % Squared Error
error_mean2=mean(error_squ2) ;  % Mean Squared Error
RMSE2(m,no_of_channels) = sqrt(error_mean2)  % Root Mean Squared Error

%[h_data2]=h_data_smooth(h_data);
figure(11*no_of_channels*m);clf; 
plot(t, 10*log10(h_data)); grid on;hold on;
%plot(t, 10*log10(h_data2)); grid on;hold on;
plot(t,10*log10(p_est_all_cls(:,m)),'color','k','LineWidth',2);
plot(t,10*log10(p_est_all),'--','color','k','LineWidth',1);
%plot(t,10*log10(h_data2),'--','color','k','LineWidth',1);
title('POWER dB')
xlabel('Time (nS)');

hold off;

end

figure(no_of_channels);clf;
plot (RMSE(:,no_of_channels))
title('Root Mean Squared Error')
xlabel('No of clusters');
set(gca,'XTick',(1:1:m));

%arrival time error;

%{
figure(2*no_of_channels);clf;
legendCell =num_of_est_clusters(:,no_of_channels)';
legendCell = cellstr(num2str(legendCell', 'No of Estimated Clusters=%-d'));
plot (abs(error_t_diffe(:,:,no_of_channels)));
xlabel('Actual Cluster ID')
ylabel('% Delay estimation difference')
legend(legendCell)
set(gca,'fontsize',12)
%}

end
 
%}