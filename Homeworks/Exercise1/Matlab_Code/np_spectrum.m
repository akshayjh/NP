function [X f] = np_spectrum(s,fs)
% INPUT
% s : Signal im Zeitbereich
% fs : Sampling-Frequenz des Signals im Zeitbereich
% 
% OUTPUT
% X : Spektrum von Zeitsignal
% f : FFT-Frequenzen

%% FFT berechnen
FFT = fft(s)/length(s);
n = length(s)/2; % H‰lfte der Koeffizienten verwerfen
X = FFT(:, 1:floor(n)+1);

%% Amplitude der Koeffizienten anpassen (auﬂer Gleichanteil)
X(2:end) = X(2:end)*2;

%% Frequenz-Vektor f berechnen
f=(0:floor(n))/n * fs/2;

%% Amplitudenspektrum plotten
figure;
plot(f,abs(X), 'LineWidth', 2), title('Spectrum of s(t), f_s = 10 kHz');
xlabel('Frequency f'), ylabel('Amplitude X');
xlim([0 10e3]);

end