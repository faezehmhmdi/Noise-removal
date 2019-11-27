clear; clc; close all;

N = 400000;
NL = 40000; %NL-> starting Noise Length

[x, fs] = audioread('src1.wav');
x = x(1:N, 1)';
x = [zeros(1, NL), x];     

%Adding uniformly distributed white noise
% noise = -0.1+rand(size(x))*(0.1-(-0.1));
y = x + 0.05*(rand(size(x)));   

figure
hold on
subplot(211);
plot(x);
title('audio sourse')
subplot(212);
plot(y);
title('noisy audio');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Computing Power spectral density of noise

frame_split = 1000; 
window_number = round(NL/frame_split);
Sz = zeros(1, frame_split/2);

freq = 0:fs/frame_split:fs/2;
freq = freq(:,2:end);

for i=0:window_number
    if i==0
        s1_fft = abs(fft(y(1:frame_split)));
        s1_fft = s1_fft(1:frame_split/2);
        s1_fft = s1_fft.*s1_fft;
        Sz1 = s1_fft/frame_split;
    else
        s1_fft = abs(fft(y(i*frame_split:(i+1)*frame_split)));
        s1_fft = s1_fft(1:frame_split/2);
        s1_fft = s1_fft.*s1_fft;
        Sz1 = s1_fft/frame_split;
    end
    Sz = Sz + Sz1;
end  

% avarage power spectral density of noise
% figure
% plot(freq,10*log10(5*Sz));
% title('average Power spectral density of noise')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Denoise algorithm

window_number = round((N+NL)/frame_split);
X = zeros(1, N);

for j=1:window_number-1
    %2.
    if j==1
        y_window = y(1:frame_split);
    else
        y_window = y((frame_split*j):(frame_split*(j+1))-1);
    end
    
    y_fft = abs(fft(y_window));
    y_fft = y_fft(1:frame_split/2);
    y_fft = y_fft.*y_fft;
    Sy = y_fft/frame_split;

    %3.
    Sx = Sy - 0.19*Sz;
    Sx(Sx < 0) = 0; 

%     if j > 100 && j<103
%         figure
%         hold on
%         plot(freq,10*log10(Sy));
%         title('Power spectral density of noisy audio')
%         xlabel('Frequency (Hz)')
%         ylabel('Power/Frequency (dB/Hz)')
%         hold off
%     end
    
    %4. 
    A1 = sqrt(Sx./Sy);
    A = [A1, fliplr(A1)];

    %5.
    Yw = fft(y_window);
    Xw = A.*Yw;
    xk = ifft(Xw);
    xk = real(xk);
    
    if j==1
        X(1:frame_split) = xk;
    else
        X((frame_split*j):(frame_split*(j+1))-1) = xk;  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Windowing denoised signal to obtain better results
window_number = round((N+NL)/frame_split);
triang_window = triang(frame_split);
triang_window = triang_window';
jump = 0.5;

x_final = zeros(1, N+NL);

for j=1:jump:window_number
    if j==1
        x_final(1:frame_split) = triang_window.*X(1:frame_split);
    else
        temp = triang_window.*X(frame_split*(j-1):frame_split*(j)-1);
        x_final(frame_split*(j-1):frame_split*(j)-1) = x_final(frame_split*(j-1):frame_split*(j)-1) + temp;
    end
end

figure
hold on
plot(y);
plot(x_final)
legend('Noisy audio','Denoised audio')
hold off

x_final = rescale(x_final, min(x), max(x));
audiowrite('denoised_audio.wav', x_final, fs);
audiowrite('noised_audio.wav', y, fs);