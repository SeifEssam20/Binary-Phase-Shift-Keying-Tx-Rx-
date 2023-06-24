clear; close all; clc;

%% Transmitter

% Define the number of bits and generate a random bit stream
no_bits = 100;
bit_stream = randi([0, 1], 1, no_bits); % 100 random bits

% Define voltage levels for line coding
V_levels = [-1 1];

% Assign voltage levels based on bit values
line_coded_stream = V_levels(bit_stream + 1);

% Initialize variables for the transmitter
Counter_1 = 1;
Counter_2 = 100;
tx_Polar_Non_Return_Zero_signal = zeros(1, 10000);

% Generate the polar non-return-to-zero (NRZ) signal
for Counter_3 = 1:no_bits
    for Counter_4 = Counter_1:Counter_2
        tx_Polar_Non_Return_Zero_signal(Counter_4) = line_coded_stream(Counter_3);
    end
    Counter_1 = Counter_1 + 100;
    Counter_2 = Counter_2 + 100;
end

% Set carrier frequency and sampling frequency
fc = 10^6;
Fs = 100 * fc;
Fs_2 = 100;
ts = 1 / Fs;

% Define the duration of each bit
Bits_duration = 0.01;

% Generate time vectors
t = 0:ts:(length(tx_Polar_Non_Return_Zero_signal) * ts) - ts; % time for the whole signal
t_2 = 0:Bits_duration:(length(tx_Polar_Non_Return_Zero_signal) * Bits_duration) - Bits_duration;

% Calculate other parameters
T = length(tx_Polar_Non_Return_Zero_signal) * ts;
T_2 = (length(tx_Polar_Non_Return_Zero_signal) * Bits_duration) ;
df = 1 / T;
df_2 = 1 /T_2;
N = ceil(T / ts);
N_2 =ceil(T_2 /Bits_duration);

% Generate the carrier signal
Carrier_Signal = cos(2 * pi * fc * t);

% Modulate the signal by multiplying with the carrier signal
modulated_signal = tx_Polar_Non_Return_Zero_signal .* Carrier_Signal;

% Check if the number of samples is even or odd for frequency calculation
if (rem(N, 2) == 0)
    freq = (-0.5 * Fs):df:(0.5 * Fs) - df; % even function
else
    freq = -(0.5 * Fs - 0.5 * df):df:(0.5 * Fs - 0.5 * df); % odd function
end

if (rem(N_2, 2) == 0)
    freq_2 = (-0.5 * Fs_2):df_2:(0.5 * Fs_2) - df_2; % even function
else
    freq_2 = -(0.5 * Fs_2 - 0.5 * df_2):df_2:(0.5 * Fs_2 - 0.5 * df_2); % odd function
end

% Plot the time domain of the polar non-return-to-zero (NRZ) signal
figure(1)
subplot(3, 1, 1)
plot(t_2, tx_Polar_Non_Return_Zero_signal);
xlabel('Time');
ylabel('Voltage');
title('Polar Non Return Zero Signal (Example)');
xlim([0 5]);
ylim([-1.5 1.5]);
grid on;

% Plot the modulated signal (BPSK)
figure(1)
subplot(3, 1, 2)
plot(t_2, modulated_signal);
xlabel('Time');
ylabel('Voltage');
title('Modulated Signal (Example)');
xlim([0 5]);
ylim([-1.5 1.5]);
grid on;

% Plot the frequency spectrum of the polar non-return-to-zero (NRZ) signal
freq_spectrum_signal = abs(fftshift(fft(tx_Polar_Non_Return_Zero_signal))).^2 / N;
figure(2)
subplot(1, 2, 1)
plot(freq_2, freq_spectrum_signal);
xlabel('Frequency');
ylabel('Magnitude squared');
title('Spectrum Signal');
grid on;

% Plot the frequency spectrum of the modulated signal
freq_spectrum_Modulated_Signal = abs(fftshift(fft(modulated_signal))).^2 / N;
figure(2)
subplot(1, 2, 2)
plot(freq, freq_spectrum_Modulated_Signal);
xlabel('Frequency');
ylabel('Magnitude squared');
title('Spectrum Modulated Signal');
grid on;


%% Receiver

% Initialize variables for the receiver
Integrated_sig =zeros(1, 100);
% Demodulate the received signal and integrate using the trapezoidal method
Recieved_sig = modulated_signal .* Carrier_Signal;

Counter_1 = 1;
Counter_2 = 100;
for Counter_3 = 1:no_bits
    for Counter_4 = Counter_1:Counter_2
        Integrated_sig(Counter_3) =Integrated_sig(Counter_3) + Recieved_sig(Counter_4) ;
    end
    Integrated_sig(Counter_3)= Integrated_sig(Counter_3) / (no_bits);
    Counter_1 = Counter_1 + 100;
    Counter_2 = Counter_2 + 100;
end

Recieved_code = Integrated_sig > 0;

% Calculate the number of errors after demodulation
Sum_errors = sum(Recieved_code ~= bit_stream);
disp('The number of errors after demodulation ='); disp(Sum_errors);

% Convert the received code back to line-coded stream
Recieved_line_coded = V_levels(Recieved_code + 1);

% Generate the received plot
Counter_1 = 1;
Counter_2 = 100;
Recieved_plot = zeros(1, 10000);

% Reconstruct the received signal
for Counter_3 = 1:no_bits
    for Counter_4 = Counter_1:Counter_2
        Recieved_plot(Counter_4) = Recieved_line_coded(Counter_3);
    end
    Counter_1 = Counter_1 + 100;
    Counter_2 = Counter_2 + 100;
end

% Plot the demodulated signal
figure(1)
subplot(3, 1, 3)
plot(t_2, Recieved_plot)
xlabel('Time');
ylabel('Voltage');
title('Demodulated Signal (Example)');
xlim([0 5]);
ylim([-1.5 1.5]);
grid on;
