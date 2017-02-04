function [s,t] = generate_signal(f_v,A_v,duration,fs)

%% Signalgeneration
% f : Frequenz-Vektor in Hz
% A : Amplituden-Vektor
% duration : Dauer des Signals in Sekunden
% fs : Abtastfrequenz in Hzf
%
sz_f = size(f_v); % Größe des Frequenz-Vektors bestimmen
t = linspace(0,duration,duration*fs); % Zeitvektor erstellen
offset = A_v(1);

%% Signal berechnen
s = zeros(1,length(t));

%%
for i = 1:sz_f(2)
    s = s + A_v(i+1)*sin(2*pi*f_v(i)*t);
end

s = s + offset;  % Gleichanteil drauf addieren

%% Signal plotten
figure;
plot(t*1e3,s, 'r');
title('Signal s(t)');
xlabel('Time t in ms'), ylabel('Amplitude A'); 
xlim([0 100]);


end