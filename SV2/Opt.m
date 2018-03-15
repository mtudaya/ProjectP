

clear all;
[h,h_ct,t,t_ct,ts]=uwb_sv_eval_ct_15_4a();
%optimization section
t_min=20;% minimum time gap between two clusters (ns)
t_min_samples=t_min/ts; % samples related to t_min
h_original=h;
h=h_original(1:t_min_samples);
tr_expected=ts;
h_len = size(h,1);
t1 = [0:(h_len-1)] * tr_expected; % 
% genrate measuremnent vectior
n=1;% number of clusters
% devide


%n_l=20; % length of the cluster
L_data= length (h);%length of the channel
n_l=floor(length (h)/n); % calculate length of the cluster
h_data=(abs(h).^2);
%h_data=h_data_mod;
h_data=h_data(1:n*n_l);
%convert h_data in to log


h_data(h_data==0) = NaN; 
h_data=log(h_data);
h_data(isnan(h_data)) = 0.00000000000000001; 

l=length(h_data)
h_data_mod=h_data;
w=ones(l,1);
for i=2:l-1
if (h_data_mod(i-1)<=h_data_mod(i)&& h_data_mod(i)>=h_data_mod(i+1))
   % h_data_mod(i)=h_data_mod(i);
    w(i)=1;
else
   % h_data_mod(i)=0;
    w(i)=0;
end
end

w=w(1:n*n_l);

% building matrix A
A = zeros((n_l*n),2*n);

beta_w=ones(n_l,1);
time_w=ones(n_l,1);
time_w = [0:(length(time_w)-1)] * (-tr_expected);
time_w=transpose(time_w);
length_diff=((n_l*n))-length(time_w);% get the length different of time_weight and column of A

for i=1:(2*n)
    if mod(i,2) == 0 
        A(:,i)=A(:,i)+[time_w;zeros(length_diff,1)];
    else
       A(:,i)= A(:,i)+[beta_w;zeros(length_diff,1)];
    end
end
cluster_arrival_intervel=n_l; % assume cluster start at 0,cluster_arrival_intervel,2*cluster_arrival_intervel

%A(:,3)=circshift(A(:,3),4)
k=0;
for i=3:2:(2*n-1)
   k=k+1;
   A(:,i)=circshift(A(:,i),k*cluster_arrival_intervel);
   A(:,i+1)=circshift(A(:,i+1),k*cluster_arrival_intervel);
   
end
y=h_data;
lambda=0.01; 
cvx_begin
            
 variable x(2*n)
            %minimize sum_square_abs(y-A*x)
  minimize( norm((w.*A * x - y), 2 ));
  %minimize( norm( A * x - h_data, 2 )+lambda* norm(x,1) );
  subject to 
        y=A*x;
        for i=1:length(x)
            if mod(i,2) == 0
            x(i)>0; 
            else
            x(i)<=0;    
            end
        end
    cvx_end
 
%exact coefficienf from solution
beta_est=zeros(n,1); % amplitude of first ray of the cluster
alpha_est=zeros(n,1);%cluster amplitude decay
for i=1:length(x)
    if mod(i,2) == 0 
        beta_est(i/2)=x(i);
    else
       alpha_est(ceil(i/2))=exp(x(i));
    end
end

%h_est=zeros(n*n_l,1);
h_est_cluster=zeros(n_l,1);
for j=1:n
for i=1:length(time_w)
   h_est_cluster(i)= alpha_est(j)*exp(beta_est(j)*time_w(i));
end
if j==1
  h_est=h_est_cluster;  
else
h_est=vertcat(h_est,h_est_cluster);
end
end
%h_est=log(h_est);
h_data=exp(h_data);



error=(h_data-h_est);    % Errors
error_squ=(h_data-h_est).^2 ;  % Squared Error
error_mean=mean(error_squ) ;  % Mean Squared Error
RMSE = sqrt(error_mean)  % Root Mean Squared Error
t=t(1:n*n_l);
figure(8); clf; plot(t, 10*log10(h_data)); grid on; hold on;
plot(t,10*log10(h_est));
title('POWER')
xlabel('Time (nS)')


