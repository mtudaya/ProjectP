
%clusers=[2,4,6,8,10,12];
clusers=[3];
for cls=1:1:length(clusers)
clearvars -except clusers cls
num_clusters_fc=clusers(cls);
[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=SV_channel_all(num_clusters_fc);

%h=10*h;
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
P_original=h_data;

%avg=mean(h_data);
min_time=2.95;%expected minimum time between clusters
%[h_data2]=h_data_smooth(h_data);
%h_data=h_data2;
h_data(h_data==0) = NaN; 
h_data(h_data<10^-15)=NaN;

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

% Find the first point where power of the channel higher than as a starting
% point
threshold=10^-10;
ix = find(h_data>threshold,1);
Y=h_data_log(ix:length(h_data));
const=1.0;

[alpha_all,beta_all,time_w_all,w_L2_all]=L2_opt2(Y,tr_expected,1,0.0001);
%[alpha_all,beta_all,time_w_all,]=L2_opt(Y,tr_expected,1);
% generate the data using alpha and beta
p_est_all=zeros(length(h_data),1);
for i=1:length(time_w_all)
   p_est_all(i)= alpha_all*exp(const*beta_all*time_w_all(i));
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
   if h_data(ray_index)> p_est_all(ray_index) || h_data(ray_index)>0.02 
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


% find the channel amplitude at the cluster starting points
%tmp_h_pwd=zeros(length(ray_index_array),1);
%tmp_h_pwd=h_data(ray_index_array);
%{
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
%}
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
  if cluster_index==1
      alpha_guess=0.0001;
  else
      alpha_guess=alpha(cluster_index-1)/10^3;
  end
      
  %[alpha(cluster_index),beta(cluster_index),time_w,w_L2]=L2_opt2(Y,tr_expected,1,alpha_guess);  
  [alpha(cluster_index),beta(cluster_index),time_w,w_L2]=L2_opt(Y,tr_expected,1);  
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

cluster_arrival_times_est=ray_index_array* (tr_expected);

p_est_all_cls(:,m) = sum(p_est_all_cluster,2); % getting the estimated power for all estimated clusters
p_est_all_cls_channel(:,m,no_of_channels) = sum(p_est_all_cluster,2);
p_est_all_channel(:,no_of_channels)=p_est_all; % estimated power of assuming one cluster
clearvars p_est_all_cluster;


%p_est_all_log=log(p_est_all); 
%p_est_all_cls_log=log(p_est_all_cls(:,m));
%no_of_clusters(m)=nnz(ray_index_array_cls_all(:,m)); % calculate number of clusters for all time gap values
%num_of_est_clusters(m,no_of_channels)=(nnz(ray_index_array_cls_all(:,m)))';% calculate number of clusters for all time gap values for all channels
%h_data(isnan(h_data)) = 0;

min_h_data=10^-15; % filter out data less than -150 dB
RMSE_weight=ones(length(h_data),1);
for weightindex=1:length(h_data)
    if h_data(weightindex)>min_h_data;
       RMSE_weight(weightindex)=1;
    else
       RMSE_weight(weightindex)=0; 
       
       
    end
end



error= RMSE_weight.*(h_data-p_est_all_cls(:,m));  % Errors
error(isnan(error))=[];
errorwwww(:,m,no_of_channels)=error;
error_squ=error.^2 ;  % Squared Error
error_squ = nonzeros(error_squ);% remove zero values of the error
errorwwww2(:,m,no_of_channels)=error_squ;
error_mean=mean(error_squ) ;  % Mean Squared Error
RMSE(m,no_of_channels) = sqrt(error_mean)  % Root Mean Squared Error

%[h_data2]=h_data_smooth(h_data);
%error2=RMSE_weight.*(h_data2-p_est_all_cls(:,m)); 
%error_squ2=error2.^2 ;  % Squared Error
%error_mean2=mean(error_squ2) ;  % Mean Squared Error
%RMSE2(m,no_of_channels) = sqrt(error_mean2)  % Root Mean Squared Error


end

for m=1:length(ray_index_array)

figure(11*no_of_channels*m);clf; 
plot(t, 10*log10(h_data)); grid on;hold on;
%plot(t, 10*log10(h_data2)); grid on;hold on;
plot(t,10*log10(p_est_all_cls(:,m)),'color','k','LineWidth',2);
plot(t,10*log10(p_est_all),'--','color','k','LineWidth',1);
%plot(t,10*log10(h_data2),'--','color','k','LineWidth',1);
title('POWER dB')
xlabel('Time (nS)');
axis([0 inf -100 10])
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

% calculate estimate cluster power





% assume number of clusters are know (
L_cls=num_clusters_fc;
Length_cluster=20;% (ns)
L_time=[0:(ceil(Length_cluster/tr_expected)-1)] * (-tr_expected);
p_est_cls_est=h_data;
h_data_new=h_data;
alpha_new=[];
beta_new=[];
alpha_new(1)=alpha(1);
beta_new(1)=beta(1);
ray_index_array_new=ray_index_array;
for cls_index=1:L_cls
  p_est_cls_est=[];  
  
  for i=1:(length(L_time))
      if i==length(h_data)
          break
      end
   p_est_cls_est(i+ray_index_array_new(cls_index)-1)= alpha_new(cls_index)*exp(beta_new(cls_index)*L_time(i));
  end
  
    p_est_cls_est=[p_est_cls_est,zeros(1,length(h_data_new)-length(p_est_cls_est))];
    p_est_cls_est_all(:,cls_index)=p_est_cls_est';
    P_cancel= h_data_new-p_est_cls_est';
    P_cancel(P_cancel==0) = NaN; 
    P_cancel(P_cancel<10^-15)=NaN;
    
    % calculate new alpha and beta
    min_P_cancel=min(P_cancel);
    P_cancel(isnan(P_cancel)) = min_P_cancel; % adjusting zero power values
    Y=log(P_cancel);


    if cls_index ~=1
        alpha_guess=alpha_new(cls_index-1)/1000;
    
        [alpha_new(cls_index+1),beta_new(cls_index+1),time_w,w_L2]=L2_opt2(Y,tr_expected,1,alpha_guess);       
    end
    
    h_data_new=[];
    h_data_new=P_cancel;
    
end
    
end


 
%filename =['X:\Google Drive\DATA\simulation\SV2\',num2str(cls,'%02d'),num2str(num_clusters_fc,'%02d'),num2str(no_of_channels,'%02d'),'data.mat']; 
%save(filename);
end

