%% Signal einlesen
%% Einlesen eins Audio Signal über den Mikrofon Eingang der Sound-Karte ###
AufnahmeZeit=2; % [s] Länge der aufgenommenen Sequenz
fs=48000;       % [Hz] Sampling Frequenz der Aufnahme
bit=16;         % [bit] Auflösung in bit
close all;
% Generierung eins Audioobjektes
recObj = audiorecorder(fs, bit, 1);
%get(recObj)
% Starten der Audioaufnahme
%disp('Start speaking.')
%pause(2); % zwei Sekunden Pause zum Luft holen :-) 
%recordblocking(recObj, AufnahmeZeit);
%disp('End of Recording.');
% Abspielen der Audio Aufnahme
%play(recObj);
% Abspeichern der Daten in ein double-precision array.
%myRecording = getaudiodata(recObj,'double');
load('myRecording.mat')
disp('Loaded soundfile from HDD.');
%save('myRecording.mat')
x = myRecording; % assign recorded signal to var 
n_ff = 44100;
%%

numChannels = 12; % test

%% Filterbank
% Greenwood Parameters (1990)
EarQ = 7.23824;
minBW = 22.8509;
order = 1;

% Berechnung der ERB & Mittenfrequenzen
% numChannels: Anzahl der Kanäle
loFreq = 100; % test
hiFreq = 8000; % test

% loFreq: Niedrigste Frequenz (zb. 100 Hz)-> Erster Filter % hiFreq: Höchste Frequenz (zb. 8000 Hz) -> Letzter Filter 
ERBlo = ((loFreq/EarQ)^order + minBW^order) ^ (1/order); 
ERBhi = ((hiFreq/EarQ)^order + minBW^order) ^ (1/order); 
overlap = (ERBhi/ERBlo)^(1/(numChannels-1));

% Equivalent Rectangular Bandwidth
ERB = ERBlo * (overlap.^(0:numChannels-1));
% Mittenfrequenzen
cf = EarQ*(((ERB.^order) - (minBW.^order)).^(1/order));

% Erzeugen der Bandpassfilter (als IIR)
ftype = 'bandpass';
f_Ny=fs/2;
N_l_Cf = (cf - ERB/2)./f_Ny;
N_u_Cf = (cf + ERB/2)./f_Ny;

for i=1:numChannels
    % Nyquist Frequenz
    % Filter Ordnung 2
    % N_l_Cf : untere Grenzfrequenz (!!Grenzfrequenz sind
    % normiert auf die Nyquist Frequenz)
    
    [b,a]=butter(2,[N_l_Cf(i) N_u_Cf(i)],ftype);
    
    % Signal Filtern
    Filterausgang(:,i)=filter(b,a,x);
    
    % Frequenzgang und Gruppenlaufzeit
    [H(:,i),w] = freqz (b, a, n_ff, fs);
    [Gr(:,i), w_gd] = grpdelay (b, a, n_ff, fs);
end
% Achtung !!!!!! Bei Filtern Höherer Ordnung (>2)
% Transfer Function Design kann instabile Filter verursachen
% (Numerische Limitation in Matlab)in diesem Fall ist ein
% Zero-Pole-Gain Design zu bevorzugen

%% Plotte alle Kanäle der fb
figure
check=0;

% for 6 channels:
if numChannels==6
for i=1:numChannels
    subplot(2,3,i), plot(linspace(0,2,length(Filterausgang)), Filterausgang(:,i))
    ylim([min(min(Filterausgang)) max(max(Filterausgang))])
    stg = ['\fontsize{18}Ch. ',num2str(i)];
    title(stg)
    if i==1
        ylabel('\fontsize{16}p in \itPa')
    end
    if i==4
        ylabel('\fontsize{16}p in \itPa')
    end
    if (mod(i,2)==1 && i >= ceil(numChannels/2))
        check=1;
    end
    if check==1
        xlabel('\fontsize{16}Time in \its')
    end
end
end


% for 12 channels:
if numChannels==12
for i=1:numChannels
    subplot(3,4,i), plot(linspace(0,2,length(Filterausgang)), Filterausgang(:,i))
    ylim([min(min(Filterausgang)) max(max(Filterausgang))])
    stg = ['\fontsize{18}Ch. ',num2str(i)];
    title(stg)
    if i==1
        ylabel('\fontsize{16}p in \itPa')
    end
    if i==5
        ylabel('\fontsize{16}p in \itPa')
    end
    if i==9
        ylabel('\fontsize{16}p in \itPa')
    end
    if (mod(i,8)==1 && i >= ceil(numChannels/2))
        check=1;
    end
    if check==1
        xlabel('\fontsize{16}Time in \its')
    end
end
end

%% Filterausgänge wieder zusammenführen und abspielen
ausgang = zeros(length(Filterausgang(:,1)), 1);
for i = 1:numChannels
    ausgang = ausgang + Filterausgang(:,i);
end
disp('Playing filtered sound\n')
pause(0.5)
sound(ausgang, fs)

%% Plots:
% Ploten des Amplitudenspektrums
% FFT des Signals
n = floor(length(x)/2);
%FFT (Normiert [0 1] bzg. Anzahl der Abtastpunkte)
FFTtmp = fft(x');
FFT = abs(FFTtmp(:,1:n+1)./n);% Real Anteil und Skalierung (Hz)
f=(0:n)/n*fs/2;
figure;
semilogx(f, FFT) % mit logarithmierter Frequenz Achse
% Darzustellender Bereich (f_min f_max FFT_min FFT_max) 
axis([50 1000 min(FFT) 0.05]); % anpassen der achsen
xlabel('\fontsize{16}Frequenz f in \itHz'); ylabel('\fontsize{16}\fontsize{16}Druck p in \itPa')
title('\fontsize{16}Einseitiges Amplitudenspektrum von y(t)');

grid on
% Ploten des Frequenzgang in KHz und dB
figure;
semilogx(w,20*log10(abs(H)), 'LineWidth', 1.5)
xlabel('\fontsize{16}Frequenz f in \itHz'); ylabel('\fontsize{16}Amplitudengang |Y(f)| in dB')
title('\fontsize{16}Einseitiges Amplitudenspektrum von y(t)');
grid on
xlim([50 10e3])
ylim([-60 0])

%% Plots Aufgabe c)
% fs: Abtastfrequenz
% Länge des Hamming Fensters in Abtastpunkten (hier 10 ms)
N_Fenster=10e-3*fs;
% Länge der Fensterüberlappung in Abtastpunkten (hier 5 ms)
Ueberlapp=5e-3*fs;
% Länge der Ausgabewerte B und f
n_out=1024;
figure
% Vor der Filterung
subplot(3,4,1:2)
plot(linspace(0,2,length(x)),x)
ylabel('p [pa]')
title('Soundfile: "Neuroprothetik"')
xlabel('Zeit in \its')
ylim([-max(abs(x)) max(abs(x))])

subplot(3,4,3:4)
plot(linspace(0,2,length(ausgang)),ausgang)
title('Gefiltert und zusammengesetzt')
xlabel('Zeit in \its')
ylim([-max(abs(ausgang)) max(abs(ausgang))])

% Vor der Filterung
[B,f,t]=specgram(x, n_out, fs, N_Fenster, Ueberlapp);
bmin= max(max(abs(B))) *10^-(60/20);
subplot(3,4,5:8)
imagesc(t, f, 20*log10( max(abs(B), bmin)/bmin ) );
ylim([0 10000])
xlabel('Zeit in \its')
ylabel('Frequenz \itf')


% Nach der Filterung
[B,f,t]=specgram(ausgang, n_out, fs, N_Fenster, Ueberlapp);
bmin= max(max(abs(B))) *10^-(60/20);
subplot(3,4,9:12)
imagesc(t, f, 20*log10( max(abs(B), bmin)/bmin ) );
ylim([0 10000])
xlabel('Zeit in \its')
ylabel('Frequenz \itf')



