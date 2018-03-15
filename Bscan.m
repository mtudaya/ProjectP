
N=60;% number of scans
M=500; %Number of samples for one scan

B=zeros(M,N);

% wall clutter from sample 50 to 100;
wall=5000*ones(1,N);

for i=50:100
    B(i,:)=1/i.*wall;
end


for i=150:170
    B(i,:)=1/i.*wall;
end


for i=250:270
    B(i,:)=1/i.*wall;
end

B(150,10)=20;
B(151,10)=20;
B(152,10)=20;
B(153,10)=20;

B(350,30)=50;
B(351,30)=50;
B(352,30)=50;
B(353,30)=50;

B(360,32)=50;
B(370,32)=50;
B(380,34)=50;
B(400,34)=50;

B(410,36)=50;
B(420,36)=50;
B(430,38)=50;
B(440,38)=50;
B(450,40)=50;
B(460,40)=50;

B(360,28)=50;
B(370,28)=50;
B(380,26)=50;
B(400,26)=50;

B(410,24)=50;
B(420,24)=50;
B(430,22)=50;
B(440,22)=50;
B(450,20)=50;
B(460,20)=50;





B(250,50)=20;
B(251,50)=20;
B(252,50)=20;
B(253,50)=20;

B(150,18)=50;
B(151,18)=50;
B(152,18)=50;
B(153,18)=50;

figure
image(B)
colormap(parula)

rank(B)

im = image(B);
axis off
colorbar
 s = size(B);
 S = s(1)*s(2);
B_vector = reshape(B,1,S);
B_vector= reshape(B',[],1);

plot (B_vector);hold on
plot(J)
J = dct2(B);

tar=J;

(abs(tar)<20);
B = A(rowInds,:);


J(abs(J) < 20) = tar(J);
J(abs(J) > 20) = 0;

K = idct2(J);
Matrix = vec2mat(K,60);
plot(Matrix)
figure
image(K)
colormap(parula)


