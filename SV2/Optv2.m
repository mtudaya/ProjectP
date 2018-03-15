
clusers=[2,4,6,8,10,12]; 
for cls=1:1:length(clusers)
clearvars -except clusers cls
num_clusters_fc=clusers(cls);
[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=uwb_sv_eval_ct_15_4a(num_clusters_fc);
%optimization section
% select all data and estimate alpha and beta asumng one cluster
min_times=[25,20,15,10,5,2,1,0.5,0.25]; % expected minimum time gap between adjucent clusters
%min_times=[25,15];
ray_index_array_cls_channel_all=[]; % estimated index of the cluster arrival for all channels with different minimum time gap (min_times) of the cluster arrivals
cluster_arrival_times= [];% estimated index of the cluster arrival with different time gap of the cluster arrivals
alpha_cls_all=[]; % estimated alpha (amplitude) for all channels 
beta_cls_all= []; % estimated beta (decay rate) for all channels 
p_est_all_channel=[];% estimated power for all channels 
p_est_all_cls_channel=[];% estimated power for all channels with different number of estinmated clusters
num_of_est_clusters=[]; % estimated number of clusters for differernt time minimum time gap of the cluster arrivals
for no_of_channels=1:num_channels
ray_index_array_cls_all=[];% estimated index of the cluster arrival  with different minimum time gap of the cluster arrivals
p_est_all_cls =[];% estimated power for all clusters 
h_data_all=(abs(h).^2);% calculation of actual channel power for all channels
h_data=(abs(h(:,no_of_channels)).^2);% calculation of actual channel power

for m=1:length(min_times)

min_time=min_times(m);%expected minimum time between clusters

h_data(h_data==0) = NaN; 
h_data_log=log(h_data);
min_h_data_log=min(h_data_log);
h_data_log(isnan(h_data)) = min_h_data_log-100; % adjusting zero power values
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
[alpha_all,beta_all,time_w,w_L2]=L2_opt(Y,tr_expected,1);
% generate the data using alpha and beta
p_est_all=zeros(n_l,1);
for i=1:length(time_w)
   p_est_all(i)= alpha_all*exp(const*beta_all*time_w(i));
end

%finding the cluster starting point (if the actual power is higher than the
%power of the one cluster model assume it is a cluster starting point

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

ray_index_array(j)=ray_index; % cluster starting point
ray_index_array_cls_all(j,m)=ray_index;
ray_index_array_cls_channel_all(j,m,no_of_channels)=ray_index;
cluster_arrival_times(j,m,no_of_channels)= ray_index* (tr_expected);
j=j+1;
ray_index=ray_index+min_sample; % shift ray index as a length of minimum cluster time gap
end
%calculate alpha and beta for clusters
alpha=zeros(length(ray_index_array),1);
beta=zeros(length(ray_index_array),1);

p_est_all_cluster=zeros(length(h_data),length(ray_index_array));
%p_est_all_cls=zeros(length(h_data),1);
%time_w_ovp = zeros(length(h_data),length(ray_index_array));

%calculate alpha and beta for clusters
for cluster_index=1:length(ray_index_array)
     
  if  cluster_index==length(ray_index_array) % check this is the last cluster, if last cluster get data from start of the cluster to end of h_data
        Y=h_data_log(ray_index_array(cluster_index):length(h_data)) ;    
  else
  
        Y=h_data_log(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1);
  end
  [alpha(cluster_index),beta(cluster_index),time_w]=L2_opt(Y,tr_expected,1);  
  alpha_cls_all(cluster_index,m,no_of_channels)=alpha(cluster_index);
  beta_cls_all(cluster_index,m,no_of_channels)=beta(cluster_index);   
 % calculate the power for a clusters   
    for k=1:length(time_w)
        p_est_all_cluster((k+ray_index_array(cluster_index)-1),cluster_index)= alpha(cluster_index)*exp(beta(cluster_index)*time_w(k));
    end

end

p_est_all_cls(:,m) = sum(p_est_all_cluster,2); % getting the estimated power for all estimated clusters
p_est_all_cls_channel(:,m,no_of_channels) = sum(p_est_all_cluster,2);
p_est_all_channel(:,no_of_channels)=p_est_all; % estimated power of assuming one cluster
clearvars p_est_all_cluster;


%p_est_all_log=log(p_est_all); 
%p_est_all_cls_log=log(p_est_all_cls(:,m));
no_of_clusters(m)=nnz(ray_index_array_cls_all(:,m)); % calculate number of clusters for all time gap values
num_of_est_clusters(m,no_of_channels)=(nnz(ray_index_array_cls_all(:,m)))';% calculate number of clusters for all time gap values for all channels
h_data(isnan(h_data)) = 0; 

min_h_data= min (h_data);
min_h_data=1.20 *min_h_data;
RMSE_weight=ones(length(h_data),1);
for weightindex=1:length(h_data)
    if h_data(weightindex)>=min_h_data;
       RMSE_weight(weightindex)=1;
    else
       RMSE_weight(weightindex)=0; 
    end
end

error=RMSE_weight.*(h_data-p_est_all_cls(:,m));    % Errors
error_squ=error.^2 ;  % Squared Error
error_mean=mean(error_squ) ;  % Mean Squared Error
RMSE(m,no_of_channels) = sqrt(error_mean)  % Root Mean Squared Error
%{
figure(11*no_of_channels*m);clf; 
plot(t, 10*log10(h_data)); grid on;hold on;
plot(t,10*log10(p_est_all_cls(:,m)),'color','r','LineWidth',1);
plot(t,10*log10(p_est_all),'--','color','k','LineWidth',1);
title('POWER dB')
xlabel('Time (nS)');
hold off;
%}
end
%{
figure(no_of_channels);clf;
plot (num_of_est_clusters(:,no_of_channels),RMSE(:,no_of_channels))
title('Root Mean Squared Error')
xlabel('No of clusters');
set(gca,'XTick',(1:1:max(num_of_est_clusters(:,no_of_channels))));
%}
%arrival time error;
[idv,idh] = size(actual_cluster_times);
[idvac,idhac,idzac] = size(cluster_arrival_times);
for mm=1:length(min_times)
    for nn=1:idv
        if(nn<=idvac)
            
            %error_t(nn,mm,no_of_channels)=(cluster_arrival_times(nn,mm,no_of_channels)-actual_cluster_times(nn,no_of_channels))/actual_cluster_times(nn,no_of_channels);
            %error_t(nn,mm,no_of_channels)=(cluster_arrival_times(nn,mm,no_of_channels)-actual_cluster_times(nn,no_of_channels));
            %error_alpha(nn,mm,no_of_channels)=(alpha_cls_all(nn,mm,no_of_channels)-Pcluster_actual(nn,no_of_channels))/Pcluster_actual(nn,no_of_channels);
            
            if actual_cluster_times(nn,no_of_channels)==0
                error_t(nn,mm,no_of_channels)=2*((cluster_arrival_times(nn,mm,no_of_channels)-actual_cluster_times(nn,no_of_channels))/(1+cluster_arrival_times(nn,mm,no_of_channels)));
                error_t_diffe(nn,mm,no_of_channels)=cluster_arrival_times(nn,mm,no_of_channels)-actual_cluster_times(nn,no_of_channels);
            else
                error_t(nn,mm,no_of_channels)=(cluster_arrival_times(nn,mm,no_of_channels)-actual_cluster_times(nn,no_of_channels))/actual_cluster_times(nn,no_of_channels);
                error_t_diffe(nn,mm,no_of_channels)=(cluster_arrival_times(nn,mm,no_of_channels)-actual_cluster_times(nn,no_of_channels));
            end
            
            if Pcluster_actual(nn,no_of_channels)==0
                error_alpha(nn,mm,no_of_channels)=((alpha_cls_all(nn,mm,no_of_channels)-Pcluster_actual(nn,no_of_channels))/(Pcluster_actual(nn,no_of_channels)+1));
            else
                error_alpha(nn,mm,no_of_channels)=(alpha_cls_all(nn,mm,no_of_channels)-Pcluster_actual(nn,no_of_channels))/Pcluster_actual(nn,no_of_channels);
            end
             
        else  % idv>idvac, actual number of clusters higher than the predicted number of clusters
            if actual_cluster_times(nn,no_of_channels)==0
                error_t(nn,mm,no_of_channels)=(actual_cluster_times(nn,no_of_channels))/(actual_cluster_times(nn,no_of_channels)+1);
                error_t_diffe(nn,mm,no_of_channels)=actual_cluster_times(nn,no_of_channels); 
            else
                error_t(nn,mm,no_of_channels)=(actual_cluster_times(nn,no_of_channels))/actual_cluster_times(nn,no_of_channels);
                error_t_diffe(nn,mm,no_of_channels)=actual_cluster_times(nn,no_of_channels);
            end
            if Pcluster_actual(nn,no_of_channels)==0
                error_alpha(nn,mm,no_of_channels)=2*Pcluster_actual(nn,no_of_channels)/(Pcluster_actual(nn,no_of_channels)+1);
            else
                error_alpha(nn,mm,no_of_channels)=Pcluster_actual(nn,no_of_channels)/Pcluster_actual(nn,no_of_channels);
            end
        
    end
end


end

error_t_squ=error_t_diffe.^2 ;  % Squared Error
error_t_mean=mean(error_t_squ);
RMSE_t(:,no_of_channels) = (sqrt(error_t_mean(:,:,no_of_channels)));  % Root Mean Squared Error
error_t_per=100*abs(error_t); % percentage error

error_alpha_squ=error_alpha.^2 ;% Squared Error
RMSE_alpha(:,no_of_channels)= sqrt(mean(error_alpha_squ(:,:,no_of_channels)));  %
error_alpha_per=100*abs(error_alpha); % percentage error
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

 
filename =['X:\Google Drive\DATA\simulation\SV2\',num2str(cls,'%02d'),num2str(num_clusters_fc,'%02d'),num2str(no_of_channels,'%02d'),'data.mat']; 
save(filename);
end
for no_of_channels=1:num_channels
plot(num_of_est_clusters(:,no_of_channels),RMSE_t(:,no_of_channels)); hold on;
end