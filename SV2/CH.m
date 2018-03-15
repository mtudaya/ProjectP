clear all;
[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=uwb_sv_eval_ct_15_4a();
%optimization section
% select all data and estimate alpha and beta asumng one cluster

for no_of_channels=1:num_channels
h_data_all=(abs(h).^2);
h_data=(abs(h(:,no_of_channels)).^2);

h_data(h_data==0) = NaN; 
h_data_log=log(h_data);
min_h_data_log=min(h_data_log);
%h_data_log(isnan(h_data)) = min_h_data_log-100; % adjusting zero power values
h_data_log(isnan(h_data)) = min_h_data_log; % adjusting zero power values
h_len = size(h_data,1);

tr_expected=ts;




figure(no_of_channels);clf; 
plot(t, 10*log10(h_data)); grid on;hold on;

colour=hsv(12);
[idv,idh] = size(actual_ray_times);
for k=1:idh
     figure(k*3); clf; 
     plot(t, 10*log10(h_data));grid on;hold on;
    for i = 1 : idv
        
      plot([actual_ray_times(i,k) actual_ray_times(i,k)], [(min(10*log10(h_data))) (max(10*log10(h_data)))],'color',colour(k,:));
       
    end
 end
%plot(t,10*log10(p_est_all_cls(:,m)),'color','r','LineWidth',2);
%plot(t,10*log10(p_est_all),'color','k','LineWidth',2);
title('POWER dB')
xlabel('Time (nS)');
hold off;

figure(no_of_channels*10); clf; 
plot(t, abs(h(:,no_of_channels))); grid on;hold on;
for idx = 1 : length(actual_cluster_times)
    plot([actual_cluster_times(idx) actual_cluster_times(idx)], [(min(abs(h(:,no_of_channels)))) (max(abs(h(:,no_of_channels))))]);
end
colour=hsv(10*num_channels);
[idv,idh] = size(actual_ray_times);
for k=1:idh
    for i = 1 : idv
     
    plot([actual_ray_times(i,k) actual_ray_times(i,k)], [(min(abs(h(:,no_of_channels)))) (max(abs(h(:,no_of_channels))))],'color',colour(k,:));
    end
 end
title('Impulse response realizations')
xlabel('Time (nS)')

figure(no_of_channels*3); clf; 
plot(t, 10*log10(h_data));grid on;hold on;
for idx = 1 : length(actual_cluster_times)
    plot([actual_cluster_times(idx) actual_cluster_times(idx)],  [(min(10*log10(h_data))) (max(10*log10(h_data)))]);
end

colour=hsv(10*num_channels);

[idv,idh] = size(actual_ray_times);
figure(no_of_channels*44); clf; 
plot(t, 10*log10(h_data));grid on;hold on;
for k=1:idh
    for i = 1 : idv
     
    plot_h=plot([actual_ray_times(i,k) actual_ray_times(i,k)],  [(min(10*log10(h_data))/k) (max(10*log10(h_data)))],'color',colour(k,:));
    end
end

title('Impulse response realizations')
xlabel('Time (nS)')

colour=hsv(10*num_channels);

[idv,idh] = size(actual_ray_times);
figure(no_of_channels*54); clf; 
plot(t, 10*log10(h_data),'LineWidth',2);grid on;hold on;
for k=1:idh
    for i = 1 : find(actual_ray_times(:,k),1,'last')
     
    plot_h=plot([actual_ray_times(i,k) actual_ray_times(i,k)],  [(min(10*log10(h_data))/k) (max(10*log10(h_data)))],'LineWidth',0.5*k,'color','k');
    end
end

title('Channel Power')
xlabel('Time (nS)')


end

