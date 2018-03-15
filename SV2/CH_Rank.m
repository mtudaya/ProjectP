
clusers=[2,4,6,8,10,12];
%clusers=[2];
e=[];
for cls=1:1:length(clusers)

num_clusters_fc=clusers(cls);
[h,h_ct,t,t_ct,ts,actual_cluster_times,Pcluster_actual,actual_ray_times,num_channels]=uwb_sv_eval_ct_15_4a(num_clusters_fc);
for no_of_channels=1:num_channels
h_data_pwd=(abs(h).^2);% calculation of actual channel power for all channels
h_data=(abs(h(:,no_of_channels)).^2);% calculation of actual channel power


r=h(:,no_of_channels);

conv_h=zeros(2*(length(r))-1,length(r));

for i=1:length(r)
    conv_h(1:length(r),i)=r;
    conv_h(:,i)=circshift(conv_h(:,i),i-1);
end
H=conv_h;
Hm=conv_h';
H_Hm=Hm*H;


end
e{cls} = eig(H_Hm);

end





