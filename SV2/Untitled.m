  Lam = 0.8;  lambda = 3;  
  lambda_mode = 0;  %0 -> Poisson process for the ray arrival times, 2 channel tap
  
  % MPC decay
  Gam =10; % higher gam low attenuation
  gamma_0 =  2; Kgamma = 0; sigma_cluster = 2.93;
  nlos = 1; 
  gamma_rise = NaN; gamma_1 = NaN; chi = NaN; % dummy in this scenario
  % Small-scale Fading
  m0 = 1; Km = 0; sigma_m0 = 1; sigma_Km = 0;
  sfading_mode = 0; m0_sp = NaN;
  % Large-scale Fading -- Shadowing
  std_shdw = 3.51;
  % Frequency Dependence
  kappa = 1.53; 
  fc = 6;   % GHz
  fs = 8;   % 2 - 10 GHz

L=1; %numberof clusters

num_channels=1;
std_L = 1/sqrt(2*Lam);      % std dev (nsec) of cluster arrival spacing
std_lam = 1/sqrt(2*lambda); % std dev (nsec) of ray arrival spacing
h_len = 1000;  % there must be a better estimate of # of paths than this
ngrow = 1000; % amount to grow data structure if more paths are needed
h = zeros(h_len,num_channels);
t = zeros(h_len,num_channels);
t0 = zeros(1,num_channels);
np = zeros(1,num_channels);
cluster_non_overlap=0; % if clusteres are not overlap 1
for k = 1:num_channels       % loop over number of channels
   
  tmp_h = zeros(size(h,1),1);
  tmp_t = zeros(size(h,1),1);
  if nlos == 1,
    Tc = (std_L*randn)^2 + (std_L*randn)^2;  % First cluster random arrival
  else
    Tc = 0;                   % First cluster arrival occurs at time 0
  end
  t0(k) = Tc;
  
  disp(L);  
  cluster_index = zeros(1,L);
  path_ix = 0;
  
  nak_m = [];

  
  for ncluster = 1:L
      
    % Determine Ray arrivals for each cluster
    Tr = 0;  % first ray arrival defined to be time 0 relative to cluster
    actual_cluster_times(ncluster,k)=Tc;
    cluster_index(ncluster) = path_ix+1 ; % remember the cluster location
    ray_ix = 0;
    gamma = Kgamma*Tc + gamma_0; % delay dependent cluster decay time
    if nlos == 2 & ncluster == 1
      gamma = gamma_1;
    end
    Mcluster = sigma_cluster*randn;
    Pcluster = 10*log10(exp(-1*Tc/Gam));%+Mcluster; % total cluster power
    Pcluster = 10^(Pcluster*0.1);
    Pcluster_actual(ncluster,k)=Pcluster; 
    
    while (Tr < 8*gamma),
      t_val = (Tc+Tr);  % time of arrival of this ray
      ray_ix=ray_ix+1;
      actual_ray_times(ray_ix,ncluster,k)=t_val;
     % if nlos == 2 & ncluster == 1
     %   h_val = Pcluster*(1-chi*exp(-Tr/gamma_rise))*exp(-Tr/gamma_1) ...
     %          *(gamma+gamma_rise)/gamma/(gamma+gamma_rise*(1-chi));
     % else
        %h_val = Pcluster/gamma*exp(-Tr/gamma);
      h_val = Pcluster*exp(-Tr/gamma)
     % end
      path_ix = path_ix + 1;  % row index of this ray
      if path_ix > h_len,
        % grow the output structures to handle more paths as needed
        tmp_h = [tmp_h; zeros(ngrow,1)];
        tmp_t = [tmp_t; zeros(ngrow,1)];
        h = [h; zeros(ngrow,num_channels)];
        t = [t; zeros(ngrow,num_channels)];
        h_len = h_len + ngrow;
      end
      tmp_h(path_ix) = h_val;
      tmp_t(path_ix) = t_val;
      if lambda_mode == 0
        Tr = Tr + (std_lam*randn)^2 + (std_lam*randn)^2;
      elseif lambda_mode == 1
        if rand < beta
          std_lam = 1/sqrt(2*lambda_1);
          Tr = Tr + (std_lam*randn)^2 + (std_lam*randn)^2;
        else
          std_lam = 1/sqrt(2*lambda_2);
          Tr = Tr + (std_lam*randn)^2 + (std_lam*randn)^2;
        end
      elseif lambda_mode == 2
        Tr = Tr + ts;
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
        %Tc = Tc + (std_L*randn)^2 + (std_L*randn)^2;
        while 1
            Tc_old=Tc;
            Tc = Tc + 10 %(std_L*randn)^2 + (std_L*randn)^2;
            if Tc > Tc_old+3;
            break
            end
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
  for path = 1:path_ix
     h_val = (gamrnd(nak_m(path), tmp_h(path)/nak_m(path))).^(1/2);
     tmp_h(path) = h_val;
  end

  np(k) = path_ix;  % number of rays (or paths) for this realization
  [sort_tmp_t,sort_ix] = sort(tmp_t(1:np(k)));    % sort in ascending time order
  t(1:np(k),k) = sort_tmp_t;
  h(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
  h_normalized(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
  % now impose a log-normal shadowing on this realization
  %fac = 10^(std_shdw*randn/20) / sqrt( h(1:np(k),k)' * h(1:np(k),k) );
  %h(1:np(k),k) = h(1:np(k),k) * fac;
 % Normalize the channel energy to 1
   h_normalized(:,k) = h_normalized(:,k)/sqrt(h_normalized(:,k)' * h_normalized(:,k) );
  % h(:,k) = h(:,k)/sqrt(h(:,k)' * h(:,k) );
   h_normalized2=h_normalized(:,k);
 actual_cluster_times2=actual_cluster_times(:,k);
 [sharedVals,idxsIntoA] = intersect(t(:,k),actual_cluster_times2);
  Pcluster_actual(:,k)=h_normalized2(idxsIntoA); 
  Pcluster_actual(:,k)= Pcluster_actual(:,k).^2;
   
end
power=abs(h).^2 
plot (t,10*log(power))
title('POWER dB')
xlabel('Time (nS)');
axis([0 inf -200 0])
return








