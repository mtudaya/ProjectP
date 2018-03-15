t=0.001:0.1:0.9;
alpha=0.5;
beta=10;
p=alpha*exp(-beta.*t);
w = 0.003 * randn(length(t),1);   
%p=p+w';
p1=real(log(p))
plot(t,p1)
Y=p1';
tr_expected=0.1;
[alpha_all,beta_all,time_w,w_L2]=L2_opt(Y,tr_expected,0);
[alpha_all2,beta_all2,time_w]=L2_opt2(Y,tr_expected);