clear all;
[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=uwb_sv_eval_ct_15_4a();
%optimization section
% select all data and estimate alpha and beta asumng one cluster

for no_of_channels=1:num_channels
h_data_all=(abs(h).^2);
h_data=(abs(h(:,no_of_channels)).^2);
%convert h_data in to log

min_times=[20,9,5,2,1,0.5];
%min_times=[229,115,60,30,20,10];
for m=1:length(min_times)
min_time=min_times(m);%minimum time between clusters
h_data(h_data==0) = NaN; 
h_data_log=log(h_data);
min_h_data_log=min(h_data_log);
h_data_log(isnan(h_data)) = min_h_data_log-100; % adjusting zero power values
h_len = size(h_data,1);
% genrate measuremnent vectior
n=1;% number of clusters
n_l=length(h_data_log); % length of the cluster
L_data= length (h_data_log);%length of the channel
tr_expected=ts;

min_sample=min_time/tr_expected;
% select the data point which are local maximum

h_data_mod_log=h_data_log;
Y=h_data_log;

[alpha_all,beta_all,time_w,w_L2]=L2_opt(Y,tr_expected,1);
% generate the data using alpha and beta
p_est_all=zeros(n_l,1);
for i=1:length(time_w)
   p_est_all(i)= alpha_all*exp(beta_all*time_w(i));
end

%finding the cluster starting point (comapring with the data of the 1
%cluster)

ray_index=1;
j=1;
ray_index_array=[];
while ray_index<length(h_data)
while true && ray_index<length(h_data)
   if h_data(ray_index)> p_est_all(ray_index);
   break ; 
   else
   ray_index=ray_index+1;
   end
end

ray_index_array(j)=ray_index
ray_index_array_cls_all(j,m)=ray_index;
cluster_arrival_times(j,m,no_of_channels)= ray_index* (tr_expected);
j=j+1;
ray_index=ray_index+min_sample; % shift ray index as a length of minimum cluster
end
%calculate alpha and beta for clusters
alpha=zeros(length(ray_index_array),1);
beta=zeros(length(ray_index_array),1);
alpha_covp=zeros(length(ray_index_array),1);
beta__covp=zeros(length(ray_index_array),1);
p_est_all_cluster=zeros(length(h_data),length(ray_index_array));
%p_est_all_cls=zeros(length(h_data),1);
time_w_ovp = zeros(length(h_data),length(ray_index_array));
for cluster_index=1:length(ray_index_array)
     
  if  cluster_index==length(ray_index_array) % check this is the last cluster, if last cluster get data from start of the cluster to end of h_data
    Y=h_data_log(ray_index_array(cluster_index):length(h_data))     
  else
  
    Y=h_data_log(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1)
  end
  [alpha(cluster_index),beta(cluster_index),time_w]=L2_opt(Y,tr_expected,1);  
    time_w_ovp(1:length(time_w),cluster_index)=time_w; 
    length_time(cluster_index)=length(time_w);
    alpha_cls_all(cluster_index,m,no_of_channels)=alpha(cluster_index);
    beta_cls_all(cluster_index,m,no_of_channels)=beta(cluster_index);   
    
    for k=1:length(time_w)
        p_est_all_cluster((k+ray_index_array(cluster_index)-1),cluster_index)= alpha(cluster_index)*exp(beta(cluster_index)*time_w(k));
    end

end
cls_overlap=2;%cluster overlap index
alpha_covp(1)=alpha(1);
beta__covp(1)=beta(1);
for cluster_index=2:length(ray_index_array)
  time_w_prv= time_w_ovp(1:length_time(cluster_index-1),cluster_index-1);
  time_w_crr= time_w_ovp(1:length_time(cluster_index),cluster_index);    
  if  cluster_index==length(ray_index_array) % check this is the last cluster, if last cluster get data from start of the cluster to end of h_data
        Y_cls_overlap=h_data_log(ray_index_array(cluster_index):length(h_data));
        time=[0:(length(Y_cls_overlap)-1)] * (-tr_expected);
        time=time+time_w_prv(length(time_w_prv));
        Y_cls_overlap2=Y_cls_overlap- 0.1*((log(alpha_covp(cluster_index-1))+ (beta__covp(cluster_index-1).*time)))';
        
         
  else
    
    Y_cls_overlap=h_data_log(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1); 
    time=[0:(length(Y_cls_overlap)-1)] * (-tr_expected);
    time=time+time_w_prv(length(time_w_prv));
    Y_cls_overlap=Y_cls_overlap-0.1*(log(alpha_covp(cluster_index-1))+ (beta__covp(cluster_index-1).*time))';
  end
    [alpha_covp(cluster_index),beta__covp(cluster_index),time_w2]=L2_opt(Y_cls_overlap,tr_expected,1);  
    
    alpha_cls_all_ovp(cluster_index,m,no_of_channels)=alpha_covp(cluster_index);
    beta_cls_all_ovp(cluster_index,m,no_of_channels)=beta__covp(cluster_index);   
    
    for k=1:length(time_w_crr)
        p_est_all_cluster_ovp((k+ray_index_array(cluster_index)-1),cluster_index)= alpha_covp(cluster_index)*exp(beta__covp(cluster_index)*time_w_crr(k));
    end   
    
    
end




p_est_all_cls(:,m) = sum(p_est_all_cluster,2);
p_est_all_cls_channel(:,m,no_of_channels) = sum(p_est_all_cluster,2);
p_est_all_channel(:,no_of_channels)=p_est_all; % estimated power of assuming one cluster
clearvars p_est_all_cluster;

p_est_all_cls_ovp(:,m) = sum(p_est_all_cluster_ovp,2);
p_est_all_cls_channel_ovp(:,m,no_of_channels) = sum(p_est_all_cluster_ovp,2);

clearvars p_est_all_cluster_ovp;


%[alpha_1,beta_1]=L2_opt(Y,tr_expected);


p_est_all_log=log(p_est_all);
p_est_all_cls_log=log(p_est_all_cls(:,m));
no_of_clusters(m)=nnz(ray_index_array_cls_all(:,m));
h_data(isnan(h_data)) = 0; 
error=(h_data-p_est_all_cls(:,m));    % Errors
error_squ=error.^2 ;  % Squared Error
error_mean=mean(error_squ) ;  % Mean Squared Error
RMSE(m,no_of_channels) = sqrt(error_mean)  % Root Mean Squared Error

figure(no_of_channels*m);clf; 
plot(t, 10*log10(h_data)); grid on;hold on;
plot(t,10*log10(p_est_all_cls(:,m)),'color','r','LineWidth',1);
%plot(t,10*log10(p_est_all_cls_ovp(:,m)),'--','color','r','LineWidth',2);
plot(t,10*log10(p_est_all),'--','color','k','LineWidth',1);
title('POWER dB')
xlabel('Time (nS)');
hold off;



end

figure(10*no_of_channels);clf;
plot (no_of_clusters,RMSE(:,no_of_channels))
title('Root Mean Squared Error')
xlabel('No of clusters');
set(gca,'XTick',(1:1:max(no_of_clusters)));

%plot(t, 10*log10(h_data)); grid on;hold on;
%plot(actual_cluster_times);


%arrival time error;
[idv,idh] = size(actual_cluster_times);
[idvac,idhac,idzac] = size(cluster_arrival_times);
for mm=1:length(min_times)
    for nn=1:idv
        if(nn<=idvac)
            error_t(nn,mm,no_of_channels)=(cluster_arrival_times(nn,mm,no_of_channels)-actual_cluster_times(nn,no_of_channels))/actual_cluster_times(nn,no_of_channels);
            error_alpha(nn,mm,no_of_channels)=(alpha_cls_all(nn,mm,no_of_channels)-Pcluster_actual(nn,no_of_channels))/Pcluster_actual(nn,no_of_channels);
        else
            error_t(nn,mm,no_of_channels)=actual_cluster_times(nn,no_of_channels)/actual_cluster_times(nn,no_of_channels);
            error_alpha(nn,mm,no_of_channels)=Pcluster_actual(nn,no_of_channels)/Pcluster_actual(nn,no_of_channels);
        
    end
end

error_t_squ=error_t.^2 ;  % Squared Error
error_alpha_squ=error_alpha.^2 ;% Squared Error
end
end
