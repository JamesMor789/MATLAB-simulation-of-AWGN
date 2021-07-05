close
clear
clc

%//////////////////////////////////////////////////////////////////////////
%Getting variable set up and generating a bit stream

N = 10;                             %number of bits
Rb = N .* 0.0001;                   %bit rate, this will lead to an Rb of 1000bps
fs=8000                             %the amount of samples to be taken
samplesPerBit = fs.*Rb;             %the number of samples to be taken per bit
y= zeros(1,N);                      %empty array to store data later
random_binary = randi([0, 1], 1,N); %create random binary data

k=1;
for i = 1:N                         %Number of bit to be sent
    for j = 1:samplesPerBit         
        y(k) = random_binary(i);    
        k=k+1;                      %this is the bit stream being sampled
    end
end


%//////////////////////////////////////////////////////////////////////////
%plotting the binary signal


bitno=1:10
figure(1)
subplot(2,1,1);
stem(bitno,random_binary)
xlabel('Bit number')
ylabel('bit high/low')
title('generated bits')

%//////////////////////////////////////////////////////////////////////////

t=linspace(0, Rb, length(y)); %time vector
fc=2000;                      %carrier frequency

figure(2)                     %plotting sampled signal in the time domain
subplot(2,1,1);
plot(t, y)
xlabel("t(s)")
ylabel("a(t)")
title("Time domain sampled signal(SNR=4)")

%/////////////////////////////////////////////////////////////////////////
%Plotting the frequency domain of the sampled signal

Num=80;
Y=fft(y);
magnitude_Y=abs(Y);    
fax_bins=[0:Num-1];    %matlab indexes starting at 1, this gets rid of that
fax_HZ=fax_bins*fs/Num; %changing the x-values froms bins to Hz
N_2 = ceil(Num/2);

figure(3)               %plotting the sampled signal(frequency domain)
plot(fax_HZ(1:N_2), magnitude_Y(1:N_2))
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('(Frequency response)Single-sided Magnitude spectrum (sampled signal(SNR=4)');

%//////////////////////////////////////////////////////////////////////////
%plotting the carrier signal and modulated signal

carrier=cos(2 .* pi .* fc * t); %carrier signal

figure(4)                       %plotting carrier signal   
plot(t, carrier)
xlabel("t(s)")
ylabel("magnitude")
title("Time domain carrier signal (SNR=4)")

modulated_signal= y .* carrier;     % modulated the signal

figure(5)                     %plotting modulated signal(time domain)
plot(t, modulated_signal)
xlabel("t(s)")
ylabel("x(t)")
title("time domain modulated signal (SNR=4)") 

Modulated_signalF=fft(modulated_signal);
magnitude_mod=abs(Modulated_signalF); 

figure(6)                     %plotting modulated signal(frequency domain)           
plot(fax_HZ(1:N_2), magnitude_mod(1:N_2))
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('(Frequency response)Single-sided Magnitude spectrum (Hz) for the modulated signal (SNR=4)');

%//////////////////////////////////////////////////////////////////////////
%creating the awgn noise and applying/plotting it to the modulated signal

transmitted_signal = awgn(modulated_signal,4,'measured') %adding noise
                                                        
                                                   
figure(7)                                  %plotting tranmitted signal                              
plot(t, transmitted_signal)
xlabel("t(s)")
ylabel("y(t)")
title("Time domain modulated signal after being transmitted (SNR=4)")

%//////////////////////////////////////////////////////////////////////////
%demodulation

recieved_signal=transmitted_signal.*carrier %demodulating the signal
                                                
figure(8)                      %plotting demodulated signal(time domain)
plot(t,recieved_signal)
xlabel("t(s)")
ylabel("y(t) post modulation")
title("Time domain demodulated signal (SNR=4)")

recieved_signalF=fft(recieved_signal);
recievedF_mod=abs(recieved_signalF);

figure(9)                    %plotting demodulated signal(frequency domain)
plot(fax_HZ(1:N_2),recievedF_mod(1:N_2))
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('(Frequency response)Single-sided Magnitude spectrum (Hz) for the recieved signal (SNR=4)')

%//////////////////////////////////////////////////////////////////////////
%filtering of the recieved signal

[b,a] = butter(5,fc/fs);              %applying thr filter to the signal
Filtered_Sig = filtfilt(b,a,recieved_signal);

figure(10)                             %plotting the filtered signal
plot(t,Filtered_Sig)
xlabel("t(s)")
ylabel("z(t)")
title("time domain filtered signal(SNR=4)") 

filtered_signalF=fft(Filtered_Sig);
filteredF_mod=abs(filtered_signalF);

figure(11)                   %plotting filtered signal(frequency domain)
plot(fax_HZ(1:N_2),filteredF_mod(1:N_2))
xlabel('Frequency(Hz)')
ylabel('Magnitude')
title('(Frequency response)Single-sided Magnitude spectrum (Hz) for the filtered signal(SNR=4)')

%//////////////////////////////////////////////////////////////////////////
%calculating a threshold value

mean_amplitude= mean(Filtered_Sig)

threshold=mean_amplitude.*0.5

%//////////////////////////////////////////////////////////////////////////
%the comparison loop using the threshold and deciding what is 1 or 0

for r = 1:80  ;                          
    
    if (Filtered_Sig(r)) >= (threshold)
        clean_signal(r) = 1;
    else
        
        clean_signal(r) = 0;
    end
end

%//////////////////////////////////////////////////////////////////////////
%the recovered bit basing its decision (0 or 1) on averaging of 8 samples.

for c=1:10
    index=(c()-1)*8+1;
    y_decod(c())=mean(clean_signal(index:index+7));
end

%//////////////////////////////////////////////////////////////////////////
%calculatinf the number of errors and BER.
  
final=round(y_decod);
numerrs = biterr(random_binary,final)  
finalBER=numerrs/N
  
figure(12)      
subplot(2,1,1);
stem(bitno,final)
xlabel('Bit number')
ylabel('bit high/low')
title('recieved bits')
  
  







