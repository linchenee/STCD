%% --- System parameter setting ---
global c0 fc lambda N K delta_f Ts Mt Mr
c0 = 3e+8;  % velocity of light
fc = 12*1e9; % carrier frequency
lambda = c0 / fc; % wavelength
N = 96; % number of subcarriers
K = 56; % number of symbols per subframe
delta_f = 240*1e3; % subcarrier spacing
T = 1 / delta_f; % effective OFDM symbol duration
Tcp = 0.5 / delta_f; % cyclic prefix duration
Ts = T + Tcp; % total OFDM symbol duration
Mt = 10; % number of transmit antennas
Mr = 10; % number of receive antennas