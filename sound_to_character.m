function sound_to_character
filename = input('Enter the file name: ')
while exist(filename) ~=2
    filename = input('File doesnot exist please Enter again: ');
end
[signal,Fs] = audioread(filename);
CharMatrix = ['1' '2' '3' 'A' '4' '5' '6' 'B' '7' '8' '9' 'C' '*' '0' '#' 'D']; % set of characters
F1 = [697 770 852 941]; % Low frequency group
F2 = [1209 1336 1477 1633];  % High frequency group
F = [F1,F2]; % All frequency array.
%%Generating Frequency Table
f  = [];
for i=1:4,
    for j=1:4,
        f = [ f ;[F1(i),F2(j)] ];
    end
end
% Fs  = 8000;       % Sampling rate = 8 kHz
t = 0:1/Fs:0.2; % length of audio signal
t_sil = 0:1/Fs:0.05; % length of zero signal
for i=1:10
    start_signal = (length(t)+length(t_sil))*(i-1)+1; % start of ith signal
    end_signal = length(t)+(length(t)+length(t_sil))*(i-1); % end of ith signal
    y(i,:)=signal(start_signal:end_signal); % segments
    
end
character = []; % to save characters and their order
L = length(t)-1; % length of the signal
nfft = length(t);
nfft2 = 2.^nextpow2(nfft); %next nearest 2.^n
for k = 1:10 % loop through the signals
    Y = fft(y(k,:)); % take fft   
    P2 = abs(Y/L); % Taking absolute value 
    P1 = P2(1:L/2+1);% to get only +ve real values by deviding axis by 2
    P1(2:end-1) = 2*P1(2:end-1);
    P1 = P1/max(P1);
    fr = Fs*(0:(L/2))/L; % define x axis in freq domain
	
    x = y(k,:);
    fx=fft(x,nfft2); %x in freq domain with signal length nfft2. gives a complex n i.e mag+phase of signal also its a mirror signal
    fx=fx(1:nfft2/2); %toget only +ve real values by deviding by 2
    Xn=fx/max(fx); %to normalize the magnitude to one

    N = 40;
    n = 0:1/Fs:N-1; % 0<n<N-1 for h(n)
    n2 = 0:29;
    L2 = length(n)-1;
    figure;
    for i=1:8
        h_n=cos(2*pi*F(i)/Fs*n); % Impulse response
        H=abs(fft(h_n)); %impulse response in freq domain
        fh1 = abs(H/n);
        fh2 = fh1(1:n/2+1);
        fh1(2:end-1) = 2*fh1(2:end-1);
        Beta=fh1/max(fh1); %to normalize the magnitude to one
        mul_a=Beta'.*P1; % filtered signal in freq domain

        h=cos(2*pi*F(i)/Fs*n2); % Impulse response
        fh=fft(h,nfft2); %impulse response in freq domain
        fh=fh(1:nfft2/2); %deviding x axis to half to get +ve real values only
        B=fh/max(fh); %to normalize the magnitude to one
        mul=B'*Xn; %filtered signal in freq domain (complex signal)
        mul1 = sum(mul,2);
        subplot(2,4,i)
        plot(abs(mul1)) %shows the output of all 8 filters for each character in seperate figure
        title(F(i)); % shows the filter freq as title
    end
    fr_comp = find(mul_a>=0.59); % find indexes greater than threshold
    fr_comp = fr(fr_comp); % get the respective frequency components
    for i=1:length(fr_comp) % to map the near frequencies values to exact value
        for j=1:length(F) % compare with each value in frequency table
            if sum(fr_comp(i) == F(j)-30:F(j)+30) % if it is in this range
                fr_comp(i) = F(j); %  map it to actuall value
            end
        end
    end
    final_freq_comp = unique(fr_comp); % remove the repeating values in array
    for i = 1:length(f) % now compare the low and high frequencies with frequency table
        if final_freq_comp == f(i,:) % if get matched
            character = [character, CharMatrix(i)]; %character will be saved
            chart = CharMatrix(i)
            soundsc(y(k,:)); % convert to sound
        end
    end
end
disp('Input Characters and their order is:');
character