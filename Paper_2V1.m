close all; clear; clc;
M = 4; % Number of antennas
B = 20e+6; % bandwidth of the receiver
Fs = 2*B; % sampling rate (2*B)
Fc = 2.4e+9; % carrier frequency
c = 3e+8; % speed of EM
lambda_em = c/Fc; % wavelength
Ts = 1/Fs; % sampling interval
L =6; % number of reflectrots
N =50; %length of the channel (-N to N)
d = lambda_em/2; % antenna spacing of the antenna array
Tau_min=0;
Tau_max=2*Ts;
ND = 8; %resolution of the dictionary
p = (Tau_max-Tau_min)/ND; %size of the dictionary
max_iteration =500;
alpha_all = zeros(L,max_iteration);
theta_all = zeros(L,max_iteration);
v_all = zeros(L,max_iteration);
Tau_estimated_all = zeros(L,max_iteration);
alpha_estimated_all = zeros(L,max_iteration);
theta_estimated_all = zeros(L,max_iteration);
v_estimated_all = zeros(L,max_iteration);
Tau_all = zeros(L,max_iteration);
error_alpha = zeros(L,max_iteration);
error_tau =zeros(L,max_iteration);
error_theta=zeros(L,max_iteration);
alpha_all_L2 = zeros(L,max_iteration);
theta_all_L2 = zeros(L,max_iteration);
v_all_L2 = zeros(L,max_iteration);
Tau_estimated_all_L2 = zeros(L,max_iteration);
alpha_estimated_all_L2 = zeros(L,max_iteration);
theta_estimated_all_L2 = zeros(L,max_iteration);
v_estimated_all_L2 = zeros(L,max_iteration);
Tau_all_L2 = zeros(L,max_iteration);
error_alpha_L2 = zeros(L,max_iteration);
error_tau_L2 =zeros(L,max_iteration);
error_theta_L2=zeros(L,max_iteration);
iteration=1;
for iteration =1:1:max_iteration
   
alpha = sort (0.1 + (1-0.1).*rand(L,1),'descend'); % r = a + (b-a).*rand(N,1)
alpha_all(:,iteration) =(alpha); 
%alpha=[0.8 0.7 0.6 0.4 0.3 0.3]; % reflection coefficients of the reflectors
%theta = [90 90 90 90 90 90]; % angle of arival of the reflectors
theta = [-30 -40 -60 -10 -20 -30];
theta_all(:,iteration) = transpose(theta);
theta_rad = deg2rad(theta);
%theta = randi([-90 90],1,L);
%v = randi([0 360],1,L);
v = [180 180 180 180 180 180]; % phase rotation due to reflection
v_all(:,iteration) = transpose(v);
v_rad=deg2rad(v);
%Tau = [(0.35*Ts) (0.55*Ts) (0.85*Ts) (2.15*Ts) (3.45*Ts) (5.65*Ts)];
%Tau = [(0.75*Ts) (1.50*Ts) (2.50*Ts) (3.50*Ts) (4.45*Ts) (5.65*Ts)];
Tau = [(0.35*Ts) (0.55*Ts) (0.75*Ts) (1.05*Ts) (1.35*Ts) (1.55*Ts)];
Tau_all(:,iteration) = transpose(Tau);
%expected delay of the reflection



 gamma = ones(M,L);
 
 for m = 1:M
        for l = 1:L
        gamma(m,l) = (2*pi/lambda_em)*d*(m-1)*(sin(theta_rad(l)));
        end
 end


size_of_samples_n = zeros(1,2*N+1);
size_of_samples_n_for_all_ant = zeros(M,2*N+1);

channel_Total = complex(size_of_samples_n_for_all_ant,0);
n = -N:1:N;
for m = 1:1:M
    
    size_of_samples_n = zeros(1,2*N+1);
    channel_RX_ANT = complex(size_of_samples_n,0);
   % channel_RX_ANT2= complex(a,0);
    %disp(m);  
    for k = 1:1:L
        %disp(k);
        %disp(channel);
        channel_RX_ANT = channel_RX_ANT + ( alpha(k)*(exp(1i*(v_rad(k)+gamma(m,k))))* (sinc(B*(n*Ts-Tau(k)-gamma(m,k)/(2*pi*Fc)))));
        %channel_RX_ANT = channel_RX_ANT + (sinc(B*((n*Ts)-T(k)-((gamma(m,k)/(2*pi*Fc))))));
        %channel_RX_ANT2 = channel_RX_ANT2 + (alpha(k)* exp(1i*(v(k)*pi/180+ gamma(m,k)))*(sinc(B*((n*Ts)-T(k)-(gamma(m,k)/(2*pi*Fc))))));
    end
    
    channel_Total(m,:)=channel_RX_ANT;
   
    
end
%plot (n,channel_RX_ANT);hold on;

% optimization section
% consider an one antenna
channel_Total_real=real(channel_Total);
channel_Total_img=imag(channel_Total);
y = channel_Total(1,:);
b= transpose(y);
% construction of dictionaty D (nxP)
D = zeros(2*N+1,ND);
% filling the D matrix

for index =-N:1:N
    
    for s = 1:1:ND
        tau_hat(s,:)=s*p;
        D((index+N+1),s)= sinc((B*((index*Ts)-(s*p))));
    end
end
DT = transpose(D);
%optimization part


%size_of_X_0ut = zeros(ND,M);
%x_BPD = complex(size_of_X_0ut,0);

%sigma = 0.0000001;                                     % Desired ||Ax - b||_2
%opts = spgSetParms('verbosity',0);

%lambda = 0.01;  
%for m = 1:1:M
% Run BPD algorithm
%[x_BPD(:,m), cost] = bpd_salsa_sparsemtx(transpose(real(channel_Total(m,:))),D, lambda, mu, Nit);
%[x_BPD(:,m)]=spg_bpdn(D,transpose((channel_Total(m,:))), sigma,opts);
%[x_BPD(:,m)]=spgl1(D,transpose((channel_Total(m,:))), 0, 1e-9, [], opts);
%[x_BPD(:,m)] = solver_L1RLS( D,transpose((channel_Total(m,:))),lambda);
%Display cost function history of BPD algorithms
%end

%for m = 1:1:M
    
%EPS=0.000001;
%mu =1;
%x0 = zeros(ND,1);
%opts = [];
%opts.maxIts     = 1000;
%opts.printEvery = 50;
%opts.tol        = 1e-8;

%lambda =0.1;
%b = transpose((channel_Total(m,:)));
%mode = 'c2c';
%AA = D; bb = b; Ahandles = linop_matrix(AA,mode);



%z0  = [];   % we don't have a good guess for the dual
%tic;
%[ x_BPD(:,m), out, optsOut ] = solver_L1RLS(Ahandles,b, lambda, x0, opts );
%time_TFOCS = toc;



%end
    for m = 1:1:M
        x0 = zeros(ND,1);
        n=1;
        y= transpose((channel_Total(m,:)));
        As=D;
        lambda =0.1;
        tic
        cvx_begin
            cvx_precision best
            variable x0(ND) complex 
            minimize lambda* norm(x0,1)+ sum_square_abs(y-As*x0)
    %minimize sum_square_abs(y-As*x0)
            subject to
            y=As*x0;
            abs( x0 ) <= 1;
        cvx_end
    [ x_BPD_CVX(:,m)]=x0;
    end 
    
  for m = 1:1:M
        x0_L2 = zeros(ND,1);
        n=1;
        y= transpose((channel_Total(m,:)));
        As=D;
        lambda =0.1;
        tic
        cvx_begin
            cvx_precision best
            variable x0_L2(ND) complex 
            minimize sum_square_abs(y-As*x0_L2)
    %minimize sum_square_abs(y-As*x0)
            subject to
            y=As*x0_L2;
            abs( x0_L2 ) <= 1;
        cvx_end
    [ x_BPD_CVX_L2(:,m)]=x0_L2;
    end 
%abs_x_BPD=abs(x_BPD);
%alpha_estimated = mean(abs_x_BPD,2);
%v_k_estimated = angle(x_BPD(:,1));
%v_k_estimated_degree =  rad2deg(v_k_estimated);
%r_2k_estimated = angle(x_BPD(:,2)) - v_k_estimated;
%r_2k_estimated_degree =  rad2deg (r_2k_estimated);
%r_3k_estimated = angle(x_BPD(:,3)) - v_k_estimated;
%r_4k_estimated = angle(x_BPD(:,4)) - v_k_estimated;
%tau_k_estimated = tau_hat -(r_2k_estimated/(2*pi*Fc));
%Sin_theta_k_estimated = r_2k_estimated *lambda_em/(2*pi*d);
%Sin_theta_k_estimated2 = gamma(2,:)*lambda_em/(2*pi*d);
%theta_k_estimated_degree = real (rad2deg((asin (Sin_theta_k_estimated))));

abs_x_BPD_CVX=abs(x_BPD_CVX);
alpha_estimated = mean(abs_x_BPD_CVX,2);
v_k_estimated = angle(x_BPD_CVX(:,1));
v_k_estimated_degree =  rad2deg(v_k_estimated);
r_2k_estimated = angle(x_BPD_CVX(:,2)) - v_k_estimated;
r_2k_estimated_degree =  rad2deg (r_2k_estimated);
r_3k_estimated = angle(x_BPD_CVX(:,3)) - v_k_estimated;
r_4k_estimated = angle(x_BPD_CVX(:,4)) - v_k_estimated;
tau_k_estimated = tau_hat -(r_2k_estimated/(2*pi*Fc));
Sin_theta_k_estimated = r_2k_estimated *lambda_em/(2*pi*d);
Sin_theta_k_estimated2 = gamma(2,:)*lambda_em/(2*pi*d);
theta_k_estimated_degree = real (rad2deg((asin (Sin_theta_k_estimated))));


%erro calcltion
[temp,originalpos] = sort( alpha_estimated, 'descend' );
positions=originalpos(1:L);
positions=sort(positions,'ascend');
alpha_selected=alpha_estimated(positions);
theta_seleted=theta_k_estimated_degree(positions);
tau_selected=tau_k_estimated(positions);
v_selected=v_k_estimated(positions);
%error_tau =
%error_theta =
error_alpha(:,iteration) = abs(alpha-alpha_selected)*100./(alpha);
error_tau(:,iteration) =abs(transpose(Tau)-tau_selected)*100./(transpose(Tau));
error_theta(:,iteration) =abs(transpose(theta)-theta_seleted)*100./(transpose(abs(theta)));
disp(iteration);

Tau_estimated_all(:,iteration) = tau_selected;
alpha_estimated_all(:,iteration) = alpha_selected;
theta_estimated_all(:,iteration) = theta_seleted;
v_estimated_all(:,iteration) = v_selected;

abs_x_BPD_CVX_L2=abs(x_BPD_CVX_L2);
alpha_estimated_L2 = mean(abs_x_BPD_CVX_L2,2);
v_k_estimated_L2 = angle(x_BPD_CVX_L2(:,1));
v_k_estimated_degree_L2 =  rad2deg(v_k_estimated_L2);
r_2k_estimated_L2 = angle(x_BPD_CVX_L2(:,2)) - v_k_estimated_L2;
r_2k_estimated_degree_L2 =  rad2deg (r_2k_estimated_L2);
r_3k_estimated_L2 = angle(x_BPD_CVX_L2(:,3)) - v_k_estimated_L2;
r_4k_estimated_L2 = angle(x_BPD_CVX_L2(:,4)) - v_k_estimated_L2;
tau_k_estimated_L2 = tau_hat -(r_2k_estimated_L2/(2*pi*Fc));
Sin_theta_k_estimated_L2 = r_2k_estimated_L2 *lambda_em/(2*pi*d);
Sin_theta_k_estimated2_L2 = gamma(2,:)*lambda_em/(2*pi*d);
theta_k_estimated_degree_L2 = real (rad2deg((asin (Sin_theta_k_estimated_L2))));


%erro calcltion
[temp,originalpos] = sort( alpha_estimated, 'descend' );
positions=originalpos(1:L);
positions=sort(positions,'ascend');
alpha_selected_L2=alpha_estimated_L2(positions);
theta_seleted_L2=theta_k_estimated_degree_L2(positions);
tau_selected_L2=tau_k_estimated_L2(positions);
v_selected_L2=v_k_estimated_L2(positions);
%error_tau =
%error_theta =
error_alpha_L2(:,iteration) = abs(alpha-alpha_selected_L2)*100./(alpha);
error_tau_L2(:,iteration) =abs(transpose(Tau)-tau_selected_L2)*100./(transpose (Tau));
error_theta_L2(:,iteration) =abs(transpose(theta)-theta_seleted_L2)*100./(transpose (abs(theta)));
disp(iteration);

Tau_estimated_all_L2(:,iteration) = tau_selected_L2;
alpha_estimated_all_L2(:,iteration) = alpha_selected_L2;
theta_estimated_all_L2(:,iteration) = theta_seleted_L2;
v_estimated_all_L2(:,iteration) = v_selected_L2;

end
difference=error_alpha-error_alpha_L2;

f1 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_alpha(1,:),'color','r','LineWidth',2);hold on;
plot (plot_x,error_alpha_L2(1,:),'color','k','LineWidth',2);
xlabel('Number of iterations')
ylabel('% Error in amplitude estimation of 1st reflector')
legend('L1','L2')
set(gca,'fontsize',16)

f2 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_alpha(2,:),'color','r','LineWidth',2);hold on;
plot (plot_x,error_alpha_L2(2,:),'color','k','LineWidth',2);
xlabel('Number of iterations')
ylabel('% Error in amplitude estimation of 2nd reflector')
legend('L1','L2')
set(gca,'fontsize',16)

f3 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_alpha(3,:),'color','r','LineWidth',2);hold on;
plot (plot_x,error_alpha_L2(3,:),'color','k','LineWidth',2);
xlabel('Number of iterations')
ylabel('% Error in amplitude estimation of 3rd reflector')
legend('L1','L2')
set(gca,'fontsize',16)

f4 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_alpha(4,:),'color','r','LineWidth',2);hold on;
plot (plot_x,error_alpha_L2(4,:),'color','k','LineWidth',2);
xlabel('Number of iterations')
ylabel('% Error in amplitude estimation of 4th reflector')
legend('L1','L2')
set(gca,'fontsize',16)

f5 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_alpha(5,:),'color','r','LineWidth',2);hold on;
plot (plot_x,error_alpha_L2(5,:),'color','k','LineWidth',2);
xlabel('Number of iterations')
ylabel('% Error in amplitude estimation of 5th reflector')
legend('L1','L2')
set(gca,'fontsize',16)

f6 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_alpha(6,:),'color','r','LineWidth',2);hold on;
plot (plot_x,error_alpha_L2(6,:),'color','k','LineWidth',2);
xlabel('Number of iterations')
ylabel('% Error in amplitude estimation of 6th reflector')
legend('L1','L2')
set(gca,'fontsize',16)


f7 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_theta(1,:),'color','r');hold on;
plot (plot_x,error_theta_L2(1,:));
f8 = figure;
plot_x=1:1:max_iteration;
plot (plot_x,error_tau(1,:),'color','r');hold on;
plot (plot_x,error_tau_L2(1,:));





      f9 = figure;
      h(1, 1) = cdfplot(error_alpha(1,:)) ;hold on; %L1 
      h(2, 1) = cdfplot(error_alpha_L2(1,:)); % L2
      xlabel('% Error in amplitude estimation 1st reflector')
     
      legend('L1','L2')
      set(gca,'fontsize',16)
      
      f10 = figure;
      h(1, 2) = cdfplot(error_alpha(2,:)) ;hold on; %%L1
      h(2, 2) = cdfplot(error_alpha_L2(2,:)); % L2
      xlabel('% Error in amplitude estimation of 2nd reflector')
      
      legend('L1','L2')
      set(gca,'fontsize',16)
      
      f11 = figure;
      h(1, 3) = cdfplot(error_alpha(3,:)) ;hold on; %L1
      h(2, 3) = cdfplot(error_alpha_L2(3,:)); % L2
      xlabel('% Error in amplitude estimation of 3rd reflector')
      
      legend('L1','L2')
      set(gca,'fontsize',16)
      
      f12 = figure;
      h(1, 4) = cdfplot(error_alpha(4,:)) ;hold on; %L1
      h(2, 4) = cdfplot(error_alpha_L2(4,:)); % L2
      xlabel('% Error in amplitude estimation of 4th reflector')
      
      legend('L1','L2')
      set(gca,'fontsize',16)
      
      f13 = figure;
      h(1, 5) = cdfplot(error_alpha(5,:)) ;hold on; %L1
      h(2, 5) = cdfplot(error_alpha_L2(5,:)); % L2
      xlabel('% Error in amplitude estimation of 5th reflector')
     
      legend('L1','L2')
      set(gca,'fontsize',16)
      
      f14 = figure;
      h(1, 6) = cdfplot(error_alpha(6,:)) ;hold on; %L1
      h(2, 6) = cdfplot(error_alpha_L2(6,:)); % L2
      xlabel('% Error in amplitude estimation of 6th reflector')
     
      legend('L1','L2')
      set(gca,'fontsize',16)
    
      
      set( h(1,:), 'LineStyle', '-', 'Color', 'r','LineWidth',2);
      set( h(2,:), 'LineStyle', '-.' , 'Color', 'k','LineWidth',2);
      
