function [alpha_est,beta_est,time_w,w_L2]=L2_opt2(y_L2,tr_expected,w_flag,alpha_guess)
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

%alpha_g=0.0001;
alpha_g=alpha_guess;
for j=1:3
%alpha_g=0.1109 % initial value
cvx_data=[];
%alpha_estimated=0.001;
for k=1:50
    %display('k' );
    %disp(k);
    fact=1/(k^2);
    cvx_begin quiet
            

        variables y b alpha_estimated
        %minimize norm (w_L2.*(y_L2 - y - (A_L2(:,2).*b)))+ 0.005*square_abs(alpha_estimated-alpha_g);
        minimize norm (w_L2.*(y_L2 - (y + A_L2(:,2).*b)))+0.05*(norm(alpha_estimated-alpha_g))
        subject to
            b>=1e-10; %power  decay rate (b) is non negative value and cannot be zero
            1e-10 <= alpha_estimated;
            alpha_estimated <= 1;
            %ln_a <= -log(alpha_g) -1-alpha_g*alpha_estimated;
            y <= log(alpha_g) + (alpha_estimated-alpha_g)/alpha_g;
    cvx_end
    cvx_data(k)=cvx_optval;
    %display('alpha_g:' );
    %disp(alpha_g);
    %display('alpha_estimated' );
    %disp(alpha_estimated);
    
    
    if square_abs(alpha_estimated-alpha_g)< 10^-10 % when norm is less than 0.001, finish this iteration
        display('error' );
        disp(square_abs(alpha_estimated-alpha_g));    
        break;
    else
        alpha_g=alpha_g+fact*(alpha_estimated-alpha_g);% if not, using the new point as the next feasible point
    
    %display('k' );
    %disp(k);
    %display('fact:' );
    %disp(fact)
   % display('error' );
   % disp(square_abs(alpha_estimated-alpha_g));
    end
end
alpha_g=alpha_estimated;
%figure(11*j);clf; 
%plot(cvx_data);

end
display('alpha_estimated' );
disp(alpha_estimated);
alpha_est=alpha_estimated;%exp(ln_a);
beta_est=b;

return; % end of program
