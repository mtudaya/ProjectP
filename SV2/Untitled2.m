
%clusers=[2,4,6,8,10,12];
clusers=[5];
for cls=1:1:length(clusers)
clearvars -except clusers cls
num_clusters_fc=clusers(cls);
[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=uwb_sv_eval_ct_15_4a(num_clusters_fc);
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
min_time=7;%expected minimum time between clusters

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

figure(11);clf; 
plot(t, 10*log10(h_data)); grid on;hold on;
[h_data2]=h_data_smooth(h_data);


end


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
 


