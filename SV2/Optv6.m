
%clusers=[2,4,6,8,10,12];
clusers=[2];
for cls=1:1:length(clusers)
clearvars -except clusers cls
num_clusters_fc=clusers(cls);
[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=uwb_sv_eval_ct_15_4a(num_clusters_fc);

h=10*h;
%optimization section
% select all data and estimate alpha and beta asumng one cluster
%min_times=[25,20,15,10,5,2,1,0.5,0.25]; % expected minimum time gap between adjucent clusters
errorwwww=[];
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
min_time=10;%expected minimum time between clusters
avg_c_length=30;% (ns)
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
tmp_h_pwd=[];
tmp_h_pwd2=[];
while ray_index<length(h_data)
    
while true && (ray_index+min_sample)<length(h_data)
   if h_data(ray_index)> p_est_all(ray_index);
        ray_index_array(j)=ray_index; % cluster starting point
   break ; 
   else
   ray_index=ray_index+1;
   end
end

%ray_index_array(j)=ray_index % cluster starting point
tmp_h_pwd(j)=h_data(ray_index);
tmp_h_pwd2(j)=p_est_all(ray_index);
%{
ray_index_array_cls_all(j,m)=ray_index;
ray_index_array_cls_channel_all(j,m,no_of_channels)=ray_index;
cluster_arrival_times(j,m,no_of_channels)= ray_index* (tr_expected);
%}
j=j+1;
ray_index=ray_index+min_sample; % shift ray index as a length of minimum cluster time gap
end

% find the channel amplitude at the cluster starting points
%tmp_h_pwd=zeros(length(ray_index_array),1);
%tmp_h_pwd=h_data(ray_index_array);

% find local maximum around the cluster staring points
for ray_index_id=1:length(ray_index_array);
temp_data=[];
if ray_index_id==length(ray_index_array)
    temp_data=h_data(ray_index_array(ray_index_id):end);
    [pks,locs] = findpeaks( temp_data);
else
    temp_data=h_data(ray_index_array(ray_index_id):ray_index_array(ray_index_id+1)-1);
    [pks,locs] = findpeaks( temp_data);
end

%[sort_temp_data, sort_temp_data_id] = sort(temp_data);
%sort_temp_data_id(end) % sort 
%ray_index_array(ray_index_id)=sort_temp_data_id(end)+ray_index_array(ray_index_id)-1;
ray_index_array(ray_index_id)=locs(1)+ray_index_array(ray_index_id)-1;
end

alpha=zeros(length(ray_index_array),1);
beta=zeros(length(ray_index_array),1);
alpha_cn=zeros(length(ray_index_array),1);
beta_cn=zeros(length(ray_index_array),1);
p_est_all_cluster=zeros(length(h_data),length(ray_index_array));
h_est_all_cluster_cn=zeros(length(h_data),length(ray_index_array));
p_est_all_cluster_cn=zeros(length(h_data),length(ray_index_array));
error_cls=[];
error_cls_all=[];
h_temp_data=zeros(length(h_data),length(ray_index_array));



%calculate alpha and beta for clusters

for cluster_index=1:length(ray_index_array)
  
    if cluster_index==1
       h_temp_data(:,cluster_index) =h_data;
        
    else
    h_temp_data(:,cluster_index)=h_data-(h_est_all_cluster_cn(:,cluster_index-1));  
    h_temp_data_negative =h_temp_data(:,cluster_index);
    h_temp_data_negative_index=find (h_temp_data_negative<0);
    h_temp_data_negative(h_temp_data_negative_index)=h_data(h_temp_data_negative_index);
    h_temp_data(:,cluster_index)=h_temp_data_negative;
    end
   
   
  if  cluster_index==length(ray_index_array) % check this is the last cluster, if last cluster get data from start of the cluster to end of h_data
        Y=h_data_log(ray_index_array(cluster_index):length(h_data)) ;
       
        Y_cn=log(h_temp_data(ray_index_array(cluster_index):length(h_temp_data),cluster_index)) ;
  elseif cluster_index==1
      
        Y=h_data_log(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1);
        Y_cn=h_data_log(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1);
       % h_temp_data =h_data(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1);
  
  else
  
        Y=h_data_log(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1);
        
        Y_cn=log(h_temp_data(ray_index_array(cluster_index):ray_index_array(cluster_index+1)-1)) ;
  end
  
  
  
  [alpha(cluster_index),beta(cluster_index),time_w,w_L2]=L2_opt(Y,tr_expected,0); 
  [alpha_cn(cluster_index),beta_cn(cluster_index),time_w_cn,w_L2_cn]=L2_opt(Y_cn,tr_expected,0); 
  
  alpha_cls_all(cluster_index,no_of_channels)=alpha(cluster_index);
  beta_cls_all(cluster_index,no_of_channels)=beta(cluster_index); 
  
  alpha_cn_cls_all(cluster_index,no_of_channels)=alpha_cn(cluster_index);
  beta_cn_cls_all(cluster_index,no_of_channels)=beta_cn(cluster_index); 
    
  
 % calculate the power for a clusters   
    for k=1:length(time_w)
        p_est_all_cluster((k+ray_index_array(cluster_index)-1),cluster_index)= alpha(cluster_index)*exp(beta(cluster_index)*time_w(k));
        p_est_all_cluster_cn((k+ray_index_array(cluster_index)-1),cluster_index)= alpha_cn(cluster_index)*exp(beta_cn(cluster_index)*time_w_cn(k));
        %calculate average error for this cluster#
        
       % error_cls(k,cluster_index)=w_L2(k)*(((h_data(k+ray_index_array(cluster_index)-1))-(alpha(cluster_index)*exp(beta(cluster_index)*time_w(k))))^2);
        error_cls(k,cluster_index)=(h_data(k+ray_index_array(cluster_index)-1)-p_est_all_cluster((k+ray_index_array(cluster_index)-1),cluster_index))^2;
        error_cls_cn(k,cluster_index)=(h_data(k+ray_index_array(cluster_index)-1)-p_est_all_cluster_cn((k+ray_index_array(cluster_index)-1),cluster_index))^2;
    end
    
    % calculate h for a clusters
    
    
    if  cluster_index==length(ray_index_array) 
        
        time_cluster=time_w_cn;
        display('time');
        disp(length(time_cluster));
        disp(cluster_index);
    else
        time_cluster=(0:ceil(avg_c_length/ts)-1)*(-1*ts);
        display('time1');
        disp(ceil(avg_c_length/ts)-1); 
        disp(cluster_index);
    end
       
    for kk=1:length(time_cluster)
        if kk+ray_index_array(cluster_index)-1==length(h)
            break
        end
        h_est_all_cluster_cn((kk+ray_index_array(cluster_index)-1),cluster_index)= alpha_cn(cluster_index)*exp(beta_cn(cluster_index)*time_cluster(kk));
              
    end   
    
    
    %calculate average error for this cluster#
  %  error_cls_sum(cluster_index,m)=sum(error_cls(:,cluster_index));
    
  % cancelation of sucssessive clusters
  
  
    
end
p_est_all_cls = sum(p_est_all_cluster,2); % getting the estimated power for all estimated clusters
p_est_all_cls_cn = sum(p_est_all_cluster_cn,2); % getting the estimated power for all estimated clusters
p_est_all_cls_channel(:,no_of_channels) = sum(p_est_all_cluster,2);
p_est_all_channel(:,no_of_channels)=p_est_all; % estimated power of assuming one cluster
%clearvars p_est_all_cluster;



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



error=(h_data-p_est_all_cls);  % Errors
errorwwww(:,no_of_channels)=error;
error_squ=error.^2 ;  % Squared Error
error_squ = nonzeros(error_squ);% remove zero values of the error
%errorwwww2(:,m,no_of_channels)=error_squ;
error_mean=mean(error_squ) ;  % Mean Squared Error
RMSE(no_of_channels) = sqrt(error_mean)  % Root Mean Squared Error

%[h_data2]=h_data_smooth(h_data);
%error2=RMSE_weight.*(h_data2-p_est_all_cls(:,m)); 
%error_squ2=error2.^2 ;  % Squared Error
%error_mean2=mean(error_squ2) ;  % Mean Squared Error
%RMSE2(m,no_of_channels) = sqrt(error_mean2)  % Root Mean Squared Error




figure(11*no_of_channels);clf; 
plot(t, 10*log10(h_data)); grid on;hold on;
%plot(t, 10*log10(h_data2)); grid on;hold on;
plot(t,10*log10(p_est_all_cls),'color','m','LineWidth',2);
plot(t,10*log10(p_est_all),'--','color','k','LineWidth',1);
plot(t,10*log10(h_est_all_cluster_cn(:,1)),'*','color','r','LineWidth',1);
plot(t,10*log10(p_est_all_cls_cn),'-.','color','r','LineWidth',1);
%plot(t,10*log10(h_data2),'--','color','k','LineWidth',1);
plot(t, 10*log10(h_data)); grid on;hold on;
plot(t,10*log10((h_temp_data(:,2))),'b--o','color','b','LineWidth',1);
title('POWER dB')
xlabel('Time (nS)');
axis([0 inf -100 10])
hold off;





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
 
%filename =['X:\Google Drive\DATA\simulation\SV2\',num2str(cls,'%02d'),num2str(num_clusters_fc,'%02d'),num2str(no_of_channels,'%02d'),'data.mat']; 
%save(filename);
end


