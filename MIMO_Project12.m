%% Initialization of Parameters
clear all
close all
c = 3*10^8;
V = 22.352;      %50 miles/hour in m/sec
BW = 30*10^3;    %Bandwidth = 30KHz
f = 800*10^6;    %Frequency  Band
lamda = c/f;
fdmax = V/lamda;
Ts = 1/BW;
alpha = 1;
k = 1:1:100000;
t = k*Ts;
N=34;
M=0.5*((0.5*N)-1);
k=1:1:100000;
SNR_db=1:20;
Mt = 2;
Mr1 = 1;
Mr2 = 2;
C_1=sqrt(2)*[exp(1i*0*(pi/2)) 0;0 exp(1i*0*(pi/2)) ];
C_2=sqrt(2)*[exp(1i*1*(pi/2)) 0;0 exp(1i*1*(pi/2)) ];
C_3=sqrt(2)*[exp(1i*2*(pi/2)) 0;0 exp(1i*2*(pi/2)) ];
C_4=sqrt(2)*[exp(1i*3*(pi/2)) 0;0 exp(1i*3*(pi/2)) ];

%% Channel Generation Using Hadamard MAtrix
H=hadamard(8);
h1_sum=zeros(1,10^5);
E1=zeros(1,8);
for q=1:1:8
    theta=(2*pi*q)/N;
    beta= pi*q/(M+1);
    for p=1:1:8
        gamma_n= (2*pi*(q+1)*p)/(M+1);
        h1 =H(q,p)*exp(1i*beta)*cos(2*pi*t'*fdmax*cos(theta) + gamma_n);
        h1_sum=h1_sum+ h1' ;
    end
    h_ch(q,:)=h1_sum;
    E1(q)=sum(power(abs(h1_sum),2))/100000;
end

for iter3=1:1:8
    h_norm(iter3,:) = h_ch(iter3,:)*(1/sqrt(E1(iter3)));
end

h1_norm = h_norm(1,:);
h2_norm = h_norm(2,:);
h3_norm = h_norm(3,:);
h4_norm = h_norm(4,:);

%% Differential MIMO System with Mr=1 & 2 fading channels
No_of_Transmissions = 10^5;
S1 = zeros(2);
I = eye(Mt);
for i = 1:1:20
    SNR=power(10,0.1*i);
    Error1 = 0;
    S0 = sqrt(Mt).*I;
    Y0=[0;0];
    for m = 1:No_of_Transmissions
        n1 = normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
        n2 = normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
        H_Ch =[h1_norm(m); h2_norm(m)];
        n = [n1 ; n2];
        
        if(m==1)
            Y0=sqrt(SNR/2).*S0*H_Ch+n;
        end
        
        S_t = (1./sqrt(2))*C_1*S0;
        Y_t = sqrt(SNR/2).*S_t*H_Ch+ n;
               
        OP_1= Y_t-(1./sqrt(2))*C_1*Y0;
        OP_2= Y_t-(1./sqrt(2))*C_2*Y0;
        OP_3= Y_t-(1./sqrt(2))*C_3*Y0;
        OP_4= Y_t-(1./sqrt(2))*C_4*Y0;
        
        Rxf1 = norm(OP_1,'fro');
        Rxf2 = norm(OP_2,'fro');
        Rxf3 = norm(OP_3,'fro');
        Rxf4 = norm(OP_4,'fro');
        
        
        if(Rxf1>Rxf2 || Rxf1>Rxf3 || Rxf1>Rxf4)
            Error1 = Error1+1;
        end
        Y0=Y_t;
        S0= S_t;
    end
    PEP1(i)= Error1/No_of_Transmissions;
end

for j=1:20
    Error2=0;
    SNR_1=power(10,0.1*j);
    S01=sqrt(2)*I;
    Y01=[0 0;0 0];
    for k = 1:No_of_Transmissions
        
        Ht_1=[h1_norm(k) h2_norm(k); h3_norm(k) h4_norm(k)];
        n1 = normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
        n2 = normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
        n3 = normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
        n4 = normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
        N1=[n1  n2; n3  n4];
        
        if(k==1)
            Y01=sqrt(SNR_1/2).*S01*Ht_1+N1;
        end
        
        St_1=(1./sqrt(2))*C_1*S01;
        Yt_1=sqrt(SNR_1/2).*St_1*Ht_1+ N1;
        
        OP_11 = Yt_1-(1./sqrt(2))*C_1*Y01;
        OP_21 = Yt_1-(1./sqrt(2))*C_2*Y01;
        OP_31 = Yt_1-(1./sqrt(2))*C_3*Y01;
        OP_41 = Yt_1-(1./sqrt(2))*C_4*Y01;
        
        Rxf11 = norm(OP_11,'fro');
        Rxf21 = norm(OP_21,'fro');
        Rxf31 = norm(OP_31,'fro');
        Rxf41 = norm(OP_41,'fro');
        
        
        if(Rxf11>Rxf21 || Rxf11>Rxf31 || Rxf11>Rxf41)
            Error2=Error2+1;
        end
        Y01=Yt_1;
        S01= St_1;
    end
    PEP2(j)=Error2/No_of_Transmissions;
end

%% Display Results
hold off;
hold on;
figure(1)
semilogy(SNR_db,PEP1,'b+',SNR_db,PEP2,'r*');
title ('Pairwise Error Probability vs SNR');
xlabel ('SNR(dB)');
ylabel ('Probability of Error');
legend('SER with Mr=1','SER with Mr=2');
grid on;
hold off;
box on





