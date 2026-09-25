%noise reduction eee243
clc;
clear;
close all;

%1.Setting Up the Necessary Parameters
Fs = 50000; % the sampling frequency in Hz
T = 1/Fs; % sampling period from Fs
t = 0:T:0.05; % time duration is from 0 to 0.05

alpha = 6; % according to group number 6
f_noise = 3870; % calculated noise frequency


%2.signal definitions

%original signal:
%x(t) = sin(2*pi*420*t)
% - cos(2*pi*190*t)
% + 0.8*sin(2*pi*250*t)*exp(-12*t)

x = sin(2*pi*420*t) ...
- cos(2*pi*190*t) ...
+ 0.8*sin(2*pi*250*t).*exp(-2*alpha*t);


%noise signal:
%n(t) = 0.5*cos(2*pi*3870*t)
% + 0.3*sin(2*pi*6500*t + 2*pi/5)

n = 0.5*cos(2*pi*f_noise*t) ...
+ 0.3*sin(2*pi*6500*t + 2*pi/5);


%noisy signal
y = x+n;


%3.FFT analysis

N = length(y);

%frequency axis
f = (0:N-1)*(Fs/N);

%positive-frequency portion
h = 1:floor(N/2);


%FFT of clean voice
X = fft(x);
Xm = abs(X)/N;


%FFT of noise
N_fft = fft(n);
Nm = abs(N_fft)/N;


%FFT of noisy signal
Y = fft(y);
Ym = abs(Y)/N;


%4.butterworth low-pass filter

filter_order = 6; %based on design
fc = 1000; %cutoff frequency in Hz

%normalized cutoff frequency
Wn = fc/(Fs/2);

%design Butterworth LPF
[b,a] = butter(filter_order,Wn,'low');

%apply zero-phase filtering
yf = filtfilt(b,a,y);

%FFT of filtered signal
Yf = fft(yf);
Yfm = abs(Yf)/N;


%5.time domain plotting

%first 0.01 seconds
tp = t <= 0.01;

figure;

subplot(4,1,1)
plot(t(tp),x(tp),'b');
title('x(t) original signal');
xlabel('time (s)');
ylabel('amplitude');
grid on;

subplot(4,1,2);
plot(t(tp),n(tp),'r');
title('n(t) noise signal');
xlabel('time (s)');
ylabel('amplitude');
grid on;

subplot(4,1,3);
plot(t(tp),y(tp),'m');
title('y(t) = x(t) + n(t) noisy signal');
xlabel('time (s)');
ylabel('amplitude');
grid on;

subplot(4,1,4);
plot(t(tp),yf(tp),'g');
title('filtered signal');
xlabel('time (s)');
ylabel('amplitude');
grid on;

sgtitle('time domain analysis');


%6.frequency domain plotting

figure;

subplot(4,1,1);
plot(f(h),2*Xm(h),'b');
title('X(f) original voice spectrum');
xlabel('frequency (Hz)');
ylabel('magnitude');
xlim([0 10000]);
grid on;

subplot(4,1,2);
plot(f(h),2*Nm(h),'r');
title('N(f) noise spectrum');
xlabel('frequency (Hz)');
ylabel('magnitude');
xlim([0 10000]);
grid on;

subplot(4,1,3);
plot(f(h),2*Ym(h),'m');
title('Y(f) noisy signal spectrum');
xlabel('frequency (Hz)');
ylabel('magnitude');
xlim([0 10000]);
grid on;

subplot(4,1,4);
plot(f(h),2*Yfm(h),'g');
title('filtered signal spectrum');
xlabel('frequency (Hz)');
ylabel('magnitude');
xlim([0 10000]);
grid on;

sgtitle('frequency domain analysis');


%7. butterworth filter frequency response

figure;

freqz(b,a,1024,Fs);

sgtitle('sixth order butterworth lowpass filter response');


%8. SNR calculation

%before filtering
SNR_before = 10*log10(mean(x.^2)/mean(n.^2));

%residual error after filtering
residual_noise = yf-x;

%after filtering
SNR_after = 10*log10(mean(x.^2)/mean(residual_noise.^2));

% improvement
SNR_cumulation = SNR_after-SNR_before;

fprintf('SNR results\n');
fprintf('SNR before filtering : %.2f dB\n',SNR_before);
fprintf('SNR after filtering : %.2f dB\n',SNR_after);
fprintf('SNR improvement : %.2f dB\n',SNR_cumulation);


%9. VERIFICATION OF THE CONTINUOUS PART OF Y(w)

%Theoretical continuous part:
%
%Yc(w) = 0.8/(2j) *
%[1/(12+j(w-500*pi)) - 1/(12+j(w+500*pi))]

f_test = [0 100 250 500 1000 2000];

fprintf('\n');

fprintf('       CONTINUOUS PART OF MANUAL Y(w) VERIFICATION\n');

fprintf(' f(Hz)       omega(rad/s)       Real{Yc}        Imag{Yc}\n');
fprintf('---------------------------------------------------------------\n');

for k = 1:length(f_test)

    w = 2*pi*f_test(k);

    Yc = (0.8/(2j)) * ...
        (1/(12 + 1j*(w-500*pi)) ...
        - 1/(12 + 1j*(w+500*pi)));

    fprintf('%5.0f       %12.4f       %12.8f    %12.8f\n', ...
        f_test(k),w,real(Yc),imag(Yc));

end



%10. INDEPENDENT NUMERICAL FOURIER TRANSFORM VERIFICATION

%Numerical CTFT of the damped 250 Hz component

Fs_verify = 200000;
dt_verify = 1/Fs_verify;

t_verify = 0:dt_verify:2;

x_damped = 0.8*sin(2*pi*250*t_verify).*exp(-12*t_verify);

fprintf('\n');
fprintf('        MANUAL FORMULA vs NUMERICAL FOURIER INTEGRAL\n');

fprintf(' f(Hz)      |Manual|       |Numerical|       Absolute Error\n');
fprintf('--------------------------------------------------------------------------\n');

for k = 1:length(f_test)

    w = 2*pi*f_test(k);

    %manual Fourier transform
    Yc_manual = (0.8/(2j)) * ...
        (1/(12 + 1j*(w-500*pi)) ...
        - 1/(12 + 1j*(w+500*pi)));

    %numerical Fourier transform
    integrand = x_damped.*exp(-1j*w*t_verify);

    Yc_numeric = trapz(t_verify,integrand);

    %absolute error
    error_value = abs(Yc_manual-Yc_numeric);

    fprintf('%5.0f       %12.8f       %12.8f       %.3e\n', ...
        f_test(k),abs(Yc_manual),abs(Yc_numeric),error_value);

end




%11. THEORETICAL FREQUENCY vs MATLAB FFT VERIFICATION

%Theoretical frequency components:
%
%190 Hz   -> -cos(2*pi*190*t)
%250 Hz   -> damped sinusoidal component
%420 Hz   -> sin(2*pi*420*t)
%3870 Hz  -> 0.5*cos(2*pi*3870*t)
%6500 Hz  -> 0.3*sin(2*pi*6500*t+2*pi/5)

%Use a longer signal for better frequency resolution

Fs_verify2 = 50000;
T_verify2 = 1/Fs_verify2;
t_verify2 = 0:T_verify2:1;

x_verify = sin(2*pi*420*t_verify2) ...
- cos(2*pi*190*t_verify2) ...
+ 0.8*sin(2*pi*250*t_verify2).*exp(-12*t_verify2);

n_verify = 0.5*cos(2*pi*3870*t_verify2) ...
+ 0.3*sin(2*pi*6500*t_verify2+2*pi/5);

y_verify = x_verify+n_verify;

N_verify = length(y_verify);

Y_verify = fft(y_verify);

f_verify = (0:N_verify-1)*(Fs_verify2/N_verify);

h_verify = 1:floor(N_verify/2);

Y_verify_mag = 2*abs(Y_verify(h_verify))/N_verify;

f_positive = f_verify(h_verify);

%theoretical frequencies
f_theory = [190 250 420 3870 6500];

fprintf('\n');

fprintf('       THEORETICAL FREQUENCY vs MATLAB FFT\n');

fprintf('Theory(Hz)    Nearest FFT(Hz)    Error(Hz)\n');


for k = 1:length(f_theory)

    [~,index] = min(abs(f_positive-f_theory(k)));

    f_found = f_positive(index);

    error_value = f_found-f_theory(k);

    fprintf('%8.0f       %12.4f       %10.4f\n', ...
        f_theory(k),f_found,error_value);

end




%12. THEORETICAL DELTA COMPONENT VERIFICATION GRAPH

%Positive-frequency delta components of Y(w):
%
%190 Hz  -> magnitude pi
%420 Hz  -> magnitude pi
%3870 Hz -> magnitude 0.5*pi
%6500 Hz -> magnitude 0.3*pi
%
%250 Hz is continuous because it is exponentially damped.

f_delta = [190 420 3870 6500];

A_delta = [pi pi 0.5*pi 0.3*pi];

figure;

stem(f_delta,A_delta,'filled');

grid on;

xlabel('frequency (Hz)');
ylabel('impulse coefficient magnitude');

title('Theoretical delta components of Y(w)');

xlim([0 7000]);


%13. THEORETICAL CONTINUOUS PART GRAPH

f_manual = linspace(0,8000,20000);

w_manual = 2*pi*f_manual;

Yc_manual_plot = (0.8/(2j)) .* ...
    (1./(12+1j*(w_manual-500*pi)) ...
    - 1./(12+1j*(w_manual+500*pi)));

figure;

plot(f_manual,abs(Yc_manual_plot),'LineWidth',1.2);

grid on;

xlabel('frequency (Hz)');
ylabel('|Yc(w)|');

title('Continuous part of theoretical Y(w)');

xlim([0 8000]);


%14. HIGH-RESOLUTION FFT VERIFICATION GRAPH

figure;

plot(f_positive,Y_verify_mag);

grid on;

xlabel('frequency (Hz)');
ylabel('magnitude');

title('High-resolution FFT of noisy signal');

xlim([0 8000]);


%15. RESIDUAL NOISE GRAPH

figure;

plot(t(tp),residual_noise(tp));

grid on;

xlabel('time (s)');
ylabel('amplitude');

title('Residual noise after filtering');


%16. RESIDUAL NOISE FREQUENCY SPECTRUM

Residual_fft = fft(residual_noise);

Residual_mag = abs(Residual_fft)/N;

figure;

plot(f(h),2*Residual_mag(h));

grid on;

xlabel('frequency (Hz)');
ylabel('magnitude');

title('Residual noise spectrum after filtering');

xlim([0 10000]);


%17. FINAL SUMMARY

fprintf('\n');

fprintf('                     PROJECT SUMMARY\n');



fprintf('noise frequency       = %d Hz\n',f_noise);
fprintf('sampling frequency    = %d Hz\n',Fs);
fprintf('cutoff frequency      = %d Hz\n',fc);

fprintf('\n');

fprintf('SNR before filtering  = %.4f dB\n',SNR_before);
fprintf('SNR after filtering   = %.4f dB\n',SNR_after);
fprintf('SNR improvement       = %.4f dB\n',SNR_cumulation);
