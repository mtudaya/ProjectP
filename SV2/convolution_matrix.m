h = [1 1 0 0 0 1 1];

r = h;
conv_h=zeros(2*(length(r))-1,length(r));

for i=1:length(r)
    conv_h(1:length(r),i)=r;
    conv_h(:,i)=circshift(conv_h(:,i),i-1);
end
H=conv_h;
Hm=conv_h';
H_Hm=Hm*H;


e= eig(H_Hm);
