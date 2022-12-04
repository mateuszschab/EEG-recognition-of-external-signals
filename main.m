%% Zadanie 1 - Widmo sygnału

clear all
close all
clc

file = ['SSVEP_96', 'SSVEP_100', 'SSVEP_88'];
load('SSVEP_88');
fs = dane_wynikowe.OpisEXP.Fs;
signal = dane_wynikowe.O1(1,:);



widmo = abs(fft(signal));
[maks, index] = max(widmo);
index_hz = (index-1) * fs/length(signal);

figure()
subplot(2,1,1)
title('Sygnał epoki I - O1')
hold on;
plot(signal)
subplot(2,1,2)
title('Widmo epoki I - O1')
hold on;
plot(widmo)

    % Wynik
 clc
 wynik = strcat("Maksymalna częstotliwość widma: ",num2str(index_hz),' Hz');
 disp('Zadanie 1:')
 disp(wynik);
 
 signal = dane_wynikowe.O1;
 fs = dane_wynikowe.OpisEXP.Fs;
for i =1:height(signal)
    widmo = abs(fft(signal(i,:)));
    [maks, index] = max(widmo);
    index_hz = (index-1) * fs/length(signal);
    recog(i,1) = index_hz;
end
    
number_correct = 0;

for i = 1:height(recog)
    if round(dane_wynikowe.klasy(i,1),2) == round(recog(i,1),2)
        number_correct = number_correct +1;
    end
end

correct_numb = (number_correct/height(recog)) * 100;

    % Wynik
    
 wynik_1 = strcat("Ilość rozpoznanych klas kanału O1: ",num2str(correct_numb),' %');
 disp('Zadanie 1:')
 disp(wynik_1);
 

%% Klasyfikator & filtr Butterwortha
clc
clear all
close all

% load('SSVEP_88');
% load('SSVEP_96');
load('SSVEP_100');

fs = dane_wynikowe.OpisEXP.Fs;
frequ = unique(dane_wynikowe.klasy);

% Wybór ilości kanałów
% ilosc_kanalow = 2;
ilosc_kanalow = 5;

for i = 1:height(frequ)
    range_freq(i,1) = frequ(i,1) - 0.5;
    range_freq(i,2) = frequ(i,1) + 0.5;
end

O1 = dane_wynikowe.O1;
O2 = dane_wynikowe.O2;
cecha = [];
warstwa = [];

for i = 1:height(O1)
    comb_O1O2(1,:) = dane_wynikowe.O1(i,:);
    comb_O1O2(2,:) = dane_wynikowe.O2(i,:);
    comb_O1O2(3,:) = dane_wynikowe.Pz(i,:);
    comb_O1O2(4,:) = dane_wynikowe.Fz(i,:);
    comb_O1O2(5,:) = dane_wynikowe.Cz(i,:);
    
     for j=1:height(range_freq)
      [a, b] = butter(4, range_freq(j, :)/(fs/2), 'bandpass');

          for k=1:ilosc_kanalow

            chosenBand = filter(a, b, comb_O1O2(k,:));
            cecha = [cecha mean(chosenBand.^2)];

          end

          warstwa(:,end+1) = cecha;
          cecha = [];

      end
    
end

Power = warstwa';
sum_power = zeros(height(Power),1);

for i =1:width(Power)
    sum_power(:,1) = sum_power(:,1) + Power(:,i);
end

class = [];

for i = 1:4:height(sum_power)
    [M,index] = max(sum_power(i:i+3));
    class(1,end+1) = frequ(index,1);
end

class = class';

number_correct = 0;

for i = 1:height(class)
    if dane_wynikowe.klasy(i,1) == class(i,1)
        number_correct = number_correct +1;
    end
end

correct_numb = (number_correct/height(class)) * 100;

    % Wynik
    
 wynik_2 = strcat("Ilość rozpoznanych klas: ",num2str(correct_numb),' %');
 disp('Zadanie 2:')
 disp(wynik_2);
 
    %% Wyniki i wnioski zadania 2
%{
    
    Sygnał 2 sekundy:
    2 kanały - 48,75%
    5 kanałów - 41,25
    
    Sygnał 3 sekundy:
    2 kanały - 68,75%
    5 kanałów - 65,00%
    
    Sygnał 4 sekundy:
    2 kanały - 65,00%
    5 kanałów - 51,25%
    
    Wnioski:
        
        Zwiększenie analizowanego czasu do trzech sekund w znacznym
    stopniu poprawiło skuteczność rozpoznania częstotliwości migotania
    diody. Jednak dalsze zwiększenie analizowanego czasu sygnału do
    czterech sekund spowodowało nie znaczy spadek skuteczności rozpoznania.
    
        Dla każdej z analizowanych długości sygnału dołączenie pozostałych
    trzech kanałów spowodowało spadek skuteczności rozpoznania
    częstotliwości migania diody.
 
    
    
%}
    
%% Zadanie 3

% unfinish

%{
clc
clear all
close all

load('SSVEP_88');
% load('SSVEP_96');
% load('SSVEP_100');

fs = dane_wynikowe.OpisEXP.Fs;
frequ = unique(dane_wynikowe.klasy);

signal = dane_wynikowe.O1(1,:);
widmo = abs(fft(signal));
[maks, index] = max(widmo);
index_hz = (index-1) * fs/length(signal);

plot(widmo);

% test = dane_wynikowe.O1(maks-5:maks+5);
N = width(dane_wynikowe.O1);
t = (0:N-1)/fs;

widmo2 = abs(fft(signal,t));
[maks, index] = max(widmo2);
index_hz = (index-1) * fs/length(signal);

plot(widmo2);


% Szukanie max. wartości widma (5)

signal = dane_wynikowe.O1(1,:);
widmo = abs(fft(signal));

for i =1:ilosc_harmonicznych
    min_dist = 1/5;
    [maks, index] = findpeaks(widmo(1:(width(widmo)/2)),'MinPeakDistance',min_dist);
    index_hz = (index-1) * fs/length(signal);
    peak_widma(1,:) = index;
    peak_widma(2,:) = index_hz;
end



N = width(signal);
t = (0:N-1)/fs;


Re=zeros(1,N);
Im=zeros(1,N);
for f=1:N
    for n=1:N
        Re(f)=Re(f)+x(n)*cos(2*pi*(f-1)*(n-1)/N);
        Im(f)=Im(f)-x(n)*sin(2*pi*(f-1)*(n-1)/N);
    end
end
Mag=(Re.^2+Im.^2).^(1/2);
Phase=atan2(Im,Re);
%}

%% Zadanie 4

clc
clear all
close all

%----------------SET VAR-------------
% load('SSVEP_88');
% load('SSVEP_96');
load('SSVEP_100');
ilosc_harmonicznych = 5;
ilosc_kanalow = 5;
%-------------------------------------

fs = dane_wynikowe.OpisEXP.Fs;
frequ = unique(dane_wynikowe.klasy);


N = width(dane_wynikowe.O1);
Ref = {};

% Do tworzenia macierzy referencji


for o = 1:height(frequ)
    r = 1;
    f = frequ(o,1);
    for i =1:ilosc_harmonicznych

        t = [0:1/fs:(N-1)/fs];
        f1 = sin(2*pi*i*f*t);
        f2 = cos(2*pi*i*f*t);

        Ref{o}(r,:) = f1;
        Ref{o}(r+1,:) = f2;
        r = r +2;
        f1 = [];
        f2 = [];
    end
end



    % Korelacja

for d = 1:height(dane_wynikowe.O1)
    
    if ilosc_kanalow == 2
        signal(1,:) = dane_wynikowe.O1(d,:);
        signal(2,:) = dane_wynikowe.O2(d,:);
    elseif ilosc_kanalow == 5
        signal(1,:) = dane_wynikowe.O1(d,:);
        signal(2,:) = dane_wynikowe.O2(d,:);
        signal(3,:) = dane_wynikowe.Pz(d,:);
        signal(4,:) = dane_wynikowe.Fz(d,:);
        signal(5,:) = dane_wynikowe.Cz(d,:);
    end
    
    
    for c = 1:width(Ref)
        [~,~,r,~,~]=canoncorr(signal',Ref{c}');
        corr(c) = r(1);
    end
    [max_val, corr_num] = max(corr);
    Correlation(d,1) = frequ(corr_num);
end
    




for c = 1:width(Correlation)
    number_correct = 0;
    for i = 1:height(Correlation)
        if dane_wynikowe.klasy(i,1) == Correlation(i,c)
            number_correct = number_correct +1;
        end
    end
    correct_numb(1,c) = (number_correct/height(Correlation)) * 100;
end


    % wynik
 disp('Zadanie 4:')
 wynik_4 = strcat("Ilość kanałów: ",num2str(ilosc_kanalow));
 disp(wynik_4);
 wynik_4 = strcat("Ilość harmonicznych: ",num2str(ilosc_harmonicznych));
 disp(wynik_4);
 wynik_4 = strcat("Ilość rozpoznanych klas dla każdego kanału: [",num2str(correct_numb),'] %');
 disp(wynik_4);

     %% Wyniki i wnioski zadania 4
%{
    
    Ilość harmonicznych: 1

Ilość kanałów: 2

2 sekundy: [68.75] %
3 sekundy: [93.75] %
4 sekundy: [91.25] %

Ilość kanałów: 5

2 sekundy: [85] %
3 sekundy: [96.25] %
4 sekundy: [100] %

Ilość harmonicznych: 2

Ilość kanałów: 2

2 sekundy:  [68.75] %
3 sekundy: [96.25] %
4 sekundy: [96.25] %

Ilość kanałów: 5

2 sekundy: [87.5] %
3 sekundy: [100] %
4 sekundy: [100] %

Ilość harmonicznych: 5

Ilość kanałów: 2

2 sekundy: [71.25] %
3 sekundy: [96.25] %
4 sekundy: [96.25] %

Ilość kanałów: 5

2 sekundy: [90] %
3 sekundy: [100] %
4 sekundy: [100] %

    Wnioski

     Zwiększenie czasu sygnału z dwóch do trzech sekund zdecydowanie
zwiększa skuteczność rozpoznawania częstotliwości. Przy 3 sekundach
rozpoznawalność jest prawie 100%, także ciężko jednoznacznie stwierdzić jak
wpływa analiza czterech sekund na poprawę skuteczności. Wartości przy
czterech sekundach generalnie pozostają bez zmian w stosunku do trzech,
jednak zdążają się pojedyncze odstępstwa, które nieznacznie wpływają na
wynik.
    
    Zwiększanie liczby harmonicznych sygnału referencyjnego poprawia
rozpoznawalność, jeżeli zwiększymy z 2 do 5 harmonicznych. Jednak
rezultaty są niewielkie np. w porównaniu do zwiększenia czasu z dwóch
do trzech sekund sygnału.
    
    Ilość kanałów, która została zwiększona z dwóch do pięciu daje podobne
rezultaty, co zwiększenie czasu analizowanego sygnału EEG.
    
%}


%% Wnioski zadań 1 - 4

%{

    Każda z wykorzystanych metod do analizy wykorzystuje inny aparat
matematyczny do osiągnięcia celu. Tym samym nie można jednoznacznie
stwierdzić, która zmiana wartości poprawia rozpoznawalność algorytmu.
Należy wziąć pod uwagę, że każde dodanie nowych wartości wymaga większej
ilości obliczeń tym samym wymagany jest dłuższy czas wykonywania operacji.
Przykładem jest tutaj funkcja korelacji, w której wprowadzenie dodatkowych
harmonicznych nieznacznie zwiększyło rozpoznawanie badanej częstotliwości
migania diody. Choć ilość danych jest niewielka to badane metody na
pojedynczych danych służą przede wszystkim przenoszenia ich na użytek w
większej skali. Dane analizowane były w trybie offline, jeżeli jednak
zaistniałaby potrzeba zaimplementowania algorytmu do analizy danych w
czasie rzeczywistym to każde dodatkowe obliczenia wprowadzałby opóźnienie w
generowanych wynikach. Dlatego biorąc pod uwagę powyższe przemyślenia,
każdy z algorytmów należy rozpatrywać pod względem skuteczności jak i czasu
obliczeń, aby znaleźć dla nich odpowiednie zastosowanie.



%}










