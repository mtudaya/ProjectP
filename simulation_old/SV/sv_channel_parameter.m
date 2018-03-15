%Function sv_channel_parameter
function [Lam,lambda,Gam,gamma,std_ln_1,std_ln_2,nlos,std_shdw] = sv_channel_parameter( channel_model )
% Return S-V model parameters for standard UWB channel models
% Lam Cluster arrival rate (clusters per nsec)
% lambda Ray arrival rate (rays per nsec)
% Gam Cluster decay factor (time constant, nsec)
% gamma Ray decay factor (time constant, nsec)
% std_ln_1 Standard deviation of log-normal variable for cluster fading
% std_ln_2 Standard deviation of log-normal variable for ray fading
% nlos Flag for non line of sight channel
% std_shdw Standard deviation of log-normal shadowing of entire impulse response
if channel_model == 1 % based on TDC measurements for LOS 0-4m
Lam = 0.0233; lambda = 2.5;
Gam = 7.1; gamma = 4.3;
std_ln_1 = 4.8 / sqrt(2);
std_ln_2 = 4.8 / sqrt(2);
nlos = 0;
std_shdw = 3;
elseif channel_model == 2 % based on TDC measurements for NLOS 0-4m
Lam = 0.4; lambda = 0.5;
Gam = 5.5; gamma = 6.7;
std_ln_1 = 4.8 / sqrt(2);
std_ln_2 = 4.8 / sqrt(2);
nlos = 1;
std_shdw = 3;
elseif channel_model == 3 % based on TDC measurements for NLOS 4-10m
Lam = 0.0667; lambda = 2.1;
Gam = 14.00; gamma = 7.9;
std_ln_1 = 4.8 / sqrt(2);
std_ln_2 = 4.8 / sqrt(2);
nlos = 1;
std_shdw = 3;
elseif channel_model == 4 % 25 nsec RMS delay spread bad multipath channel
Lam = 0.0667; lambda = 2.1;
Gam = 24; gamma = 12;
std_ln_1 = 4.8 / sqrt(2);
std_ln_2 = 4.8 / sqrt(2);
nlos = 1;
std_shdw = 3;
end
return
 end