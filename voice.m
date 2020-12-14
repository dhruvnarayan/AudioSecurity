clc;
clear;
close all;

% Uncomment Code to get audio input -> Either 1. Record audio file or 2. Load audio file
% 1. Record your voice for 5 seconds.
% recObj = audiorecorder(8000,8,2);   % 8000 Hz - Sample Rate Fs, 8-bit, 2 channels
% %disp('Start speaking.')
% recordblocking(recObj,5);     
% %disp('End of Recording.')
% 
% %Play back the input audio file.
% %play(recObj);                           
% 
% %Store data in double-precision array.
% N = getaudiodata(recObj);
% 
% %Saving N to a MAT file - 
% save('example.mat', 'N');

% 2. Load audio file
% load('example.mat');

% Plotting the important figures
subplot(311);
plot(N);
title('Time domain characteristic of Input Plain Audio');
xlabel('Time(ms)');
ylabel('Amplitude');
axis([0 40000 -1 1]);

% Start Timer
tic;

% L is length of the dual channel matrix (rows of Input_Data)
L = size(N,1);

% Begin Hashing of Plain Audio file

Opt.method = 'SHA-256';
Opt.format = 'hex';
Hash1 = DataHash(N,Opt);

%Hashing of Audio File with one bit change
M = N;
if N(1,1)~=0
    M(1,1) = - N(1,1);
else
    M(1,1) = 0.5;
end

Opt.method = 'SHA-256';
Opt.format = 'hex';
Hash2 = DataHash(M,Opt);

% Cryptosystem -1

%Calculating initial conditions
x_key = 1.1 * 10^-14;
y_key = 1.2 * 10^-14;
z_key = 1.3 * 10^-14;
x_initial = x_key + hex2dec(Hash1(1:8))/(L*(10^14));
y_initial = y_key + hex2dec(Hash1(9:16))/(L*(10^14));
z_initial = z_key + hex2dec(Hash1(17:24))/(L*(10^14));
transform_len = 2 * mod(hex2dec(Hash1(25:28)),6) + 3;
st = 500 + mod(hex2dec(Hash1(29:32)),L)/1000;

% initial conditions for runge-kutta
x(1) = x_initial;
y(1) = y_initial;
z(1) = z_initial;
w(1) = 4.4;
h = 0.002;    % step size
a = 10;
b = 8/3;
c = 28;
r = -0.5;  %-1.52 < r < -0.06

% Calculation loop for Runge-Kutta
for i = 1:(L/2)
    
    %Runge-Kutta solution for x - trajectory
    k11 = (a * (y(i)-x(i)))+w(i);
    k12 = (a * (y(i)-(x(i)+ h*k11/2))) + w(i);
    k13 = (a * (y(i)-(x(i)+ h*k12/2))) + w(i);
    k14 = (a * (y(i)-(x(i)+ h*k13))) + w(i);
    x(i+1) = x(i)+((h/6)*(k11+k12+k13+k14));
    
    %Runge-Kutta solution for y - trajectory
    k21 = c*x(i)-y(i)-x(i)*z(i);
    k22 = c*x(i)-(y(i)+k21*h/2)-x(i)*z(i);
    k23 = c*x(i)-(y(i)+k22*h/2)-x(i)*z(i);
    k24 = c*x(i)-(y(i)+k23*h/2)-x(i)*z(i);
    y(i+1) = y(i)+((h/6)*(k21+k22+k23+k24));
    
    %Runge-Kutta solution for z- trajectory
    k31 = x(i)*y(i) - b*z(i);
    k32 = x(i)*y(i) - b*(z(i)+k31*h/2);
    k33 = x(i)*y(i) - b*(z(i)+k32*h/2);
    k34 = x(i)*y(i) - b*(z(i)+k33*h);
    z(i+1) = z(i)+((h/6)*(k31+k32+k33+k34));
    
    %Runge-Kutta solution for w
    k41 = -y(i)*z(i) + r*w(i);
    k42 = -y(i)*z(i) + r*(w(i)+k41*h/2);
    k43 = -y(i)*z(i) + r*(w(i)+k42*h/2);
    k44 = -y(i)*z(i) + r*(w(i)+k43*h);
    w(i+1) = w(i)+((h/6)*(k11+k12+k13+k14));
    
end

% Cryptosystem -2

%Calculating initial conditions
xx_key = 2.1 * 10^-14;
yy_key = 2.2 * 10^-14;
zz_key = 2.3 * 10^-14;
xx_initial = x_key + hex2dec(Hash2(1:8))/(L*(10^14));
yy_initial = y_key + hex2dec(Hash2(9:16))/(L*(10^14));
zz_initial = z_key + hex2dec(Hash2(17:24))/(L*(10^14));
transform_len_2 = 2 * mod(hex2dec(Hash2(25:28)),6) + 3;
st_2 = 500 + mod(hex2dec(Hash2(29:32)),L)/1000;

% initial conditions for runge-kutta
xx(1) = x_initial;
yy(1) = y_initial;
zz(1) = z_initial;
ww(1) = 4.4;
%h = 0.002;    % step size
% a = 10;
% b = 8/3;
% c = 28;
% r = -0.5;  %-1.52 < r < -0.06


% Calculation loop for Runge-Kutta
for i = ((L/2)+1):L
    
    %Runge-Kutta solution for xx - trajectory
    k11 = (a * (y(i)-x(i)))+w(i);
    k12 = (a * (y(i)-(x(i)+ h*k11/2))) + w(i);
    k13 = (a * (y(i)-(x(i)+ h*k12/2))) + w(i);
    k14 = (a * (y(i)-(x(i)+ h*k13))) + w(i);
    x(i+1) = x(i)+((h/6)*(k11+k12+k13+k14));
    
    %Runge-Kutta solution for yy - trajectory
    k21 = c*x(i)-y(i)-x(i)*z(i);
    k22 = c*x(i)-(y(i)+k21*h/2)-x(i)*z(i);
    k23 = c*x(i)-(y(i)+k22*h/2)-x(i)*z(i);
    k24 = c*x(i)-(y(i)+k23*h/2)-x(i)*z(i);
    y(i+1) = y(i)+((h/6)*(k21+k22+k23+k24));
    
    %Runge-Kutta solution for zz- trajectory
    k31 = x(i)*y(i) - b*z(i);
    k32 = x(i)*y(i) - b*(z(i)+k31*h/2);
    k33 = x(i)*y(i) - b*(z(i)+k32*h/2);
    k34 = x(i)*y(i) - b*(z(i)+k33*h);
    z(i+1) = z(i)+((h/6)*(k31+k32+k33+k34));
    
    %Runge-Kutta solution for ww
    k41 = -y(i)*z(i) + r*w(i);
    k42 = -y(i)*z(i) + r*(w(i)+k41*h/2);
    k43 = -y(i)*z(i) + r*(w(i)+k42*h/2);
    k44 = -y(i)*z(i) + r*(w(i)+k43*h);
    w(i+1) = w(i)+((h/6)*(k11+k12+k13+k14));
    
end

% Line numbers of dual-channel matrix N in ascending order
H1 = 1:L;
H2 = 1:L;

% Permuting rows of H1 and H2
for i = L:-1:2
    H11(i) = H1(mod(floor(abs(x(i))*L),i)+1);
    H22(i) = H2(mod(floor(abs(z(i))*L),i)+1);
end  

% Getting sequences of Yc and Zc
for i=1:L
    Yc(i) = mod(10*y(i),1);
    Zc(i) = mod(z(i),1);
end

%Generating Confusion Matrixx - Nh
Nh = zeros(L,2);
for i=1:L
    Nh(i,1) = N(H11(i)+1,1);
    Nh(i,2) = N(H22(i)+1,1);
end

% Diffusion process
Nc = zeros(L,2);
for i=1:L
    Nc(i,1)=(Yc(i) * (1 + (Yc(i) * Nh(i,1)) + 1 - Nh(i,1)))/ b;
    Nc(i,2)=(Zc(i) * (1 + (Zc(i) * Nh(i,2)) + 1 - Nh(i,2)))/ b;
end
 

% Plotting the important figures
subplot(312);
plot(Nc);
title('Time domain characteristic of Ciphered Audio');
xlabel('Time(ms)');
ylabel('Amplitude');
axis([0 40000 0 1]);

%ENCRYPTION COMPLETE - 'Nc' contains encrypted signal


%DECRYPTION PROCESS BEGIN

%The line numbers H1 and H2 are same for Nc as N (declared before)

%xx,y,z vectors generated by Runge-Kutta method remain the same

%Permuting rows of H1 and H2 to get H11 and H22 (already declared above)

%Yc and Zc remain the same

%Inversely diffuse the elements in Nc to get dec_Nh
dec_Nh = zeros(L,2);
for i = 1:L
    dec_Nh(i,1) = (b* Nc(i,1) - 2*Yc(i)) / (Yc(i)*(Yc(i)-1));
    dec_Nh(i,2) = (b* Nc(i,2) - 2*Zc(i)) / (Zc(i)*(Zc(i)-1));
end

%Inversely permute two rows of dec_Nh to get deciphered audio data matrixx
dec_N = zeros(L,2);
for i= 1:L
    dec_N(H11(i)+1,1) = dec_Nh(i,1);
    dec_N(H22(i)+1,2) = dec_Nh(i,2);
end

%DECRYPTION END - dec_N containts the deciphered audio matrixx

% Plotting the important figures
subplot(313);
plot(dec_N);
title('Time domain characteristic of Deciphered Audio');
xlabel('Time(ms)');
ylabel('Amplitude');
axis([0 40000 -1 1]);

% Stop Timer
toc;

% Determining PSNR between N and Nc

N_min = abs(min(min(N)))+N;
Nc_min = abs(min(min(Nc)))+Nc;

temp = max(max([N_min;Nc_min]));

N_new = round((N_min/temp)*255);
Nc_new = round((Nc_min/temp)*255);

[psnr,mse,maxxerr] = psnr_mse_maxerr11(N_new,Nc_new)


% Determining correlation between channels of N

N_left = N(:,1);
N_right = N(:,2);
r_1 = corrcoef(N_left,N_right)

% Determining correlation between channels of Nc

Nc_left = Nc(:,1);
Nc_right = Nc(:,2);
r_2 = corrcoef(Nc_left,Nc_right)

% Determining correlation between channels of N and Nc

N_left = N(:,1);
Nc_left = Nc(:,1);
r_3 = corrcoef(N_left,Nc_left)

N_right = N(:,2);
Nc_right = Nc(:,2);
r_4 = corrcoef(N_right,Nc_right)
