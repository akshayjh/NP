%% Signal einlesen
%% Einlesen eins Audio Signal über den Mikrofon Eingang der Sound-Karte ###
clear all, clc
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

numChannels = 16;
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


%% Berechnung der Einhüllenden
% [s] Hilbert Trans. über das Ausgangssignal
HT = hilbert(Filterausgang); % [pa] Komplexwertige Amplituden
Envelope = abs(HT);
% [pa] Betrag der Amplituden
% filtern der Filterausgänge mit Tiefpass @30 Hz Eckfrequenz
ftype = 'low';
for i=1:numChannels
    [b,a]=butter(2, 30/fs/2,ftype);
    Envelope_filtered(:,i) = filter(b,a,Envelope(:,i));
    
end

%% Original Sound:
% Vor der Filterung
figure
plot(linspace(0,2,length(x)),x, 'r')
ylabel('\fontsize{18}p in \itPa')
title('\fontsize{18}Soundfile: "Neuroprothetik"')
xlabel('\fontsize{18}Zeit in \its')
ylim([-max(abs(x))*1.2 max(abs(x))*1.2])


%% 1a ) plot Filterausgänge

if numChannels == 8
    figure
    for i = 1:numChannels
        subplot(2,4,i), plot(linspace(0,2,length(x)),Envelope_filtered(:,i), '-r', 'LineWidth', 1.5)
        if i==1 || i == 5
            ylabel('\fontsize{16}p in \itPa');
        end
        if i>=5
           xlabel('\fontsize{16}Zeit t in \its'); 
        end
        stg=['\fontsize{16}Channel Nr.', num2str(i)];
        title(stg);
        %xlim([0 fs/10])
        grid on
    end
end


if numChannels == 16
    figure
    for i = 1:numChannels
        subplot(4,4,i), plot(linspace(0,2,length(x)),Envelope_filtered(:,i), '-r', 'LineWidth', 1.5)
        if i == 1 || i == 5 || i == 9 || i == 13
            ylabel('\fontsize{16}p in \itPa');
        end
        if i>=13
           xlabel('\fontsize{16}Zeit t in \its'); 
        end
        stg=['\fontsize{16}Channel Nr.', num2str(i)];
        title(stg);
        %xlim([0 fs/10])
        grid on
    end
end

%% 1b) plot 
% Ploten des Amplitudenspektrums
figure;
subplot(1,2,1)
imagesc(Envelope_filtered');
title('\fontsize{16}Envelope')
ylabel('\fontsize{16}Ch. Nr.')
xlabel('\fontsize{16}Zeit t in \its');

% gemittelt über die Zeit
subplot(1,2,2)
Envelope_mean = mean(Envelope_filtered);
bar(Envelope_mean);
title('\fontsize{16}Mean of each channel')
ylabel('\fontsize{16}p in \itPa')
xlabel('\fontsize{16}Ch. Nr.');


%% Dynamikkompression
SPL_Schwelle=30; % Hoerschwelle in [dB]
SPL_Dynamic=50; % Dynamikbereich in [dB]
I_THR=250e-3;     % [A] Threshold
I_MCL=10*I_THR;    % [A] Most Confortabel Level


%---------------------------------------------------
% A_log=500; % Kompensation MedEl : 250 500 750 1000 
%---------------------------------------------------
c=1000; % Compression rate

%[pa] Untere _Dynamik Schwelle
    p_min=10^( (SPL_Schwelle) / 20 ) * 20e-6; %[pa]
%[pa] Obere _Dynamik Schwelle
    p_max=10^( (SPL_Schwelle+SPL_Dynamic) / 20 ) * 20e-6;%[pa]
    
% Map Rule
% Akustisches Signal normieren
p_normed = (Envelope_filtered - p_min) / (p_max - p_min); % [0 1]
% "Clipping" der Werte außerhalb des Mapping Bereichs vermeiden
p_normed(p_normed>1) = 1;
    
if c == 0
% c = 0 bedeutet lineares Mapping (siehe Handbuch MedEL)
% [pa] -> [0 1] -> [A]
    Envelops_Dy_Comp = (I_MCL - I_THR) * p_normed + I_THR;
else
    % [pa] -> [0 1] -> [A]
    Envelops_Dy_Comp = (I_MCL - I_THR) * log10(1 + c * p_normed) / log10(c + 1) + I_THR;
end
% Unterschwellige Werte auf 0 setzen
Envelops_Dy_Comp (Envelops_Dy_Comp < I_THR*2) = 0;
Envelops_Dy_Comp = real(Envelops_Dy_Comp);

%% plot 
figure

% for 8-channel case
if numChannels == 8
    for i=1:numChannels
        subplot(2,4,i), plot(Envelops_Dy_Comp(:,i), '-r', 'LineWidth', 1.5)
        x = ['\fontsize{16}Ch. Nr.', num2str(i)];
        title(x)
        if i == 1 || i == 5
            ylabel('\fontsize{16}Stimulus in mA');
        end
        if i >= 5
            xlabel('\fontsize{16}Zeit in \its')
        end
    end
    grid on
end

% for 16-channel case
if numChannels == 16
    for i=1:numChannels
        subplot(4,4,i), plot(Envelops_Dy_Comp(:,i), '-r', 'LineWidth', 1.5)
        x = ['\fontsize{16}Ch. Nr.', num2str(i)];
        title(x)
        if i == 1 || i == 5 || i == 9 || i == 13
            ylabel('\fontsize{16}Stimulus in mA');
        end
        if i >= 13
            xlabel('\fontsize{16}Zeit in \its')
        end
    end
    grid on
end


%% 2b)
% specgram
figure
subplot(1,2,1)
imagesc(Envelops_Dy_Comp');
title('\fontsize{16}Freq. Kanäle in Falschfarben')
xlabel('\fontsize{16}Zeit'), ylabel('\fontsize{16}Ch. Nr')

% gemittelt über die Zeit
subplot(1,2,2)
Envelope_dy_mean = mean(Envelops_Dy_Comp);
bar(Envelope_dy_mean);
title('\fontsize{16}Mean of each channel')
ylabel('\fontsize{16}p in \itPa')
xlabel('\fontsize{16}Ch. Nr.');

%% 2c)
figure
tit = [num2str(SPL_Dynamic),'\fontsize{16} [dB] Dynamic'];
subplot(1,2,1), plot(Envelope_filtered(:,3), Envelops_Dy_Comp(:,3), '-g', 'LineWidth', 2), grid on;
title(tit)
ylabel('\fontsize{16}I in \itmA')
xlabel('\fontsize{16}p in \itPa')
grid on
grid minor

subplot(1,2,2), semilogx(Envelope_filtered(:,8), Envelops_Dy_Comp(:,8), '-g', 'LineWidth', 2), grid on;
title(tit)
ylabel('\fontsize{16}I in \itmA')
xlabel('\fontsize{16}p in \itPa')

grid on 
grid minor
xlim([10e-5 1])



