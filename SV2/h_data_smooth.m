function [h_data_sm]=h_data_smooth(h_data)

y_smooth=log(h_data);
%y_smooth=h_data;
ln_y_smooth=length(y_smooth);
%plot((y))
ln_y_smooth=ln_y_smooth-1;
Pi=zeros(ln_y_smooth-1,ln_y_smooth);


d1 = -1*ones(1,ln_y_smooth);
d0 = ones(1,ln_y_smooth);
Ad0=diag(d0);
Ad0=[Ad0 zeros(size(Ad0,1),1)];
Ad1=diag(d1,1) ;
Ad1(end,:) = [];
Pi = Ad0+Ad1 ;
Pi2=Pi;
Pi2(end,:) = [];
Pi2=Pi2(:,1:end-1);
MaxIt=5;
W_smooth=ones((ln_y_smooth-1),1);    % initial weights
W_smooth_m=diag(W_smooth);

epsilen=0.05;
for i=1:MaxIt
    W_new2(:,i)=W_smooth;
    cvx_begin quiet
            
        variable h_data_sm(length(y_smooth))
            %minimize sum_square_abs(y-A*x)
  
        minimize( norm( y_smooth-h_data_sm, 2 ));
        subject to 
            norm(W_smooth_m*Pi2*Pi*h_data_sm,1)<=20;
       
     cvx_end
     
     W_coefficient=abs(h_data_sm(1:end-2))+epsilen;
     W_smooth=1./W_coefficient;
     %big = y_smooth <-70;
     W_smooth(big) = false;
     W_smooth_m=diag(W_smooth);
    
end
 figure(110);clf; 
 h_PI=Pi2*Pi*h_data_sm;
 h_data_sm=exp(h_data_sm);
 plot((y_smooth));hold on;
 plot(log(h_data_sm),'color','k','LineWidth',2);hold on;
% plot(h_PI ,'LineWidth',1);hold on;

  
   return 