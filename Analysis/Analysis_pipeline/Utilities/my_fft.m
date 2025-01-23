%% spectrum function
function [P1, f] = my_fft(X, Fs, L)
%     Fs = 1/dt *1000.;
%     freq = 0:Fs/Ttot:Fs/2;
%     %df = Fs/Ttot;
% 
%     zf0 = fft(zz - mean(zz));
%     zf = abs(zf0(1:length(zf0)/2)) / (Ttot/(Fs/1e3));
%     zf(2:end) = 2*zf(2:end);

dt = 1/Fs;             % Sampling period       
t = (0:L-1)*dt;        % Time vector
f = Fs*(0:(L/2))/L;    % frequency

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

end