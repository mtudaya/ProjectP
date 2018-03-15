N=length(h_data);
Pi1=zeros(N-2,N-1);
Pi2=zeros(N-1,N);

d1 = 1:5 ;
d0 = 1:6 ;
d2 = 1:5 ;
A = diag(d2,-1)+diag(d0)+diag(d1,1) ;

lambd=0.1;
cvx_begin quiet
            
 variable P_hat(N)
            
  minimize( norm( h_data - P_hat, 2 )+lambda* norm(Pi1*Pi2*P_hat,1) );
  subject to 
        Pi1*Pi2*P_hat<15;
        
    cvx_end