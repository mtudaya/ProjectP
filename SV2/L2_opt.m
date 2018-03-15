function [alpha_est,beta_est,time_w,w_L2]=L2_opt(y_L2,tr_expected,w_flag)
L_data=length(y_L2);
w=ones(L_data,1);


if w_flag==1 && L_data>1
    for i=1:L_data
        if i==1 
           
            if (y_L2(i)>=y_L2(i+1))
            % h_data_mod(i)=h_data_mod(i);
                w(i)=1;
                
            else
            % h_data_mod(i)=0;
            w(i)=0;
            end
       
        elseif i==L_data
            if (y_L2(i-1)<=y_L2(i))
         % h_data_mod(i)=h_data_mod(i);
                w(i)=1;
            else
         % h_data_mod(i)=0;
                w(i)=0;
            end
        
        else
     
            if (y_L2(i-1)<=y_L2(i)&& y_L2(i)>=y_L2(i+1))
        % h_data_mod(i)=h_data_mod(i);
                w(i)=1;
            else
             % h_data_mod(i)=0;
                w(i)=0;
            end
        end
        if y_L2(i)< -65
            w(i)=0;
        end
       
    end
       
elseif w_flag==0 
    w=ones(L_data,1); 
elseif  L_data==1
    w=ones(L_data,1); 
else
    disp(w_flag);
    disp(L_data);
    error('w_flag numner wrong!!')
    
end
% building matrix A
A = zeros(L_data,2);

beta_w=ones(L_data,1);
time_w=ones(L_data,1);
time_w = [0:(length(time_w)-1)] * (-tr_expected);
time_w=transpose(time_w);
A(:,2)=time_w;
A(:,1)=beta_w;
w_L2=w;
A_L2=A;

cvx_begin quiet
            

    variables ln_a b
    minimize(norm (w_L2.*(y_L2 - (ln_a + A_L2(:,2).*b))))
    subject to
        b>=1e-5; %power  decay rate (b) is non negative value and cannot be zero
        ln_a<=0
cvx_end
alpha_est=exp(ln_a);
beta_est=b;

return; % end of program
