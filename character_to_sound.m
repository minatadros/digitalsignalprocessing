function character_to_sound

CharMatrix = ['1' '2' '3' 'A' '4' '5' '6' 'B' '7' '8' '9' 'C' '*' '0' '#' 'D']; % set of characters
F1 = [697 770 852 941]; % Low frequency group
F2 = [1209 1336 1477 1633];  % High frequency group
%%Generating Frequency Table
f  = [];
for i=1:4,
    for j=1:4,
        f = [ f ;[F1(i),F2(j)] ];
    end
end
%% a)  First it should ask user to enter ten  characters
x=input('Please enter 10 Characters, Alphabets should be capital: ');  % Enter matrix of 10 characters from given Character set e.g. ['1' '6' 'D' '8' '*'.....upto 10 chars
while length(x)~=10 || sum(ismember(x, CharMatrix))~=10;
    if length(x)~=10
        x=input('Input Characters are not exact 10 Characters. Please Try Again, Alphabets should be capital'); % warn the user to enter characters from the valid set
    elseif sum(ismember(x, CharMatrix))~=10;
        x=input('Input Characters are not valid. Please Try Again, Alphabets should be capital'); % warn the user to enter characters from the valid set
    end
end
% x = '1234567890';
%% b)  For  each character,  your code should produce  summation  of  two  sinusoids
Fs  = 8000;       % Sampling rate = 8 kHz
t = 0:1/Fs:0.2; % time for signal
t_sil = 0:1/Fs:0.05;    % time for silence signal
sil = zeros(1,length(t_sil)); % zero signal
signal = []; % to concatinate the all signals
for i=1:10
    k = find(CharMatrix == x(i)); % find the character in CharMatrix
    freq_table(i,:) = f(k,:); % pick their corresponding frequencies
    y(i,:) = 0.5*cos(2*pi*freq_table(i,1)*t) + 0.5*cos(2*pi*freq_table(i,2)*t); % sum up the signals
    signal = [signal y(i,:) sil]; % conactenate them
end
soundsc(signal); % convert to sound
filename = 'character_to_sound_data.wav'; % write to .wav file
audiowrite(filename,signal,Fs)