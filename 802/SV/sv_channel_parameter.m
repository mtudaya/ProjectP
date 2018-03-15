function [Lam_cluster,lambda_ray,Gam_cluster,gamma_ray,std_cluster_fading_db,std_ray_fading_db,nlos,std_shadow] = sv_channel_parameter( channel_num )
%Lam_cluster Cluster arrival rate (clusters per nsec)
%lambda_ray  Ray arrival rate (rays per nsec)
%Gam_cluster Cluster decay factor (time constant, nsec)
%gamma_ray  Ray decay factor (time constant, nsec)
%std_cluster_fading_db Standard deviation of log-normal variable for cluster fading 
%std_ray_fading_db Standard deviation of log-normal variable for ray fading
%std_shadow Standard deviation of log-normal shadowing of entire impulse response


if channel_num == 1 % based on TDC measurements for LOS 0-4m
Lam_cluster = 0.0233; lambda_ray = 2.5;
Gam_cluster = 7.1; gamma_ray = 4.3;
std_cluster_fading_db = 4.8 / sqrt(2);
std_ray_fading_db = 10.6 / sqrt(2);
nlos = 0;
std_shadow = 3;
elseif channel_num == 2 % based on TDC measurements for NLOS 0-4m
Lam_cluster = 0.4; lambda_ray = 0.5;
Gam_cluster = 5.5; gamma_ray = 6.7;
std_cluster_fading_db = 4.8 / sqrt(2);
std_ray_fading_db = 4.8 / sqrt(2);
nlos = 1;
std_shadow = 3;
elseif channel_num == 3 % based on TDC measurements for NLOS 4-10m
Lam_cluster = 0.0667; lambda_ray = 2.1;
Gam_cluster = 14.00; gamma_ray = 7.9;
std_cluster_fading_db = 4.8 / sqrt(2);
std_ray_fading_db = 4.8 / sqrt(2);
nlos = 1;
std_shadow = 3;
elseif channel_num == 4 % 25 nsec RMS delay spread bad multipath channel
Lam_cluster = 0.0667; lambda_ray = 2.1;
Gam_cluster = 24; gamma_ray = 12;
std_cluster_fading_db = 4.8 / sqrt(2);
std_ray_fading_db = 4.8 / sqrt(2);
nlos = 1;
std_shadow = 3;
end
return