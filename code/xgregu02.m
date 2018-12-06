%%% First task %%%
[s, Fs] = audioread('xgregu02.wav');
period = 16;  % this is one period
%fileID = fopen('tmp_xgregu02.txt', 'w');
k = 1;
i = 9;

fprintf("1. Úloha\n");
fprintf("Vzorkovacia frekvencia je: %f [Hz].\n", Fs);
fprintf("Dĺžka signálu vo vzorkoch: %f, v sekundách: %f [s].\n", length(s), length(s)/Fs);

fprintf("\n2. Úloha\n");
%%% Second task %%%
% Decoding the input into signals
for i=1:period:size(s)
  if(s(i+8) > 0)
    decoded(k) = 1;
  elseif(s(i+8) < 0)
    decoded(k) = 0;
  end
  k = k + 1;
end

%%% The Symbols %%%
x_axis_dec = 0.0005:0.001:0.020;
y_decoded = decoded(1:20);
%%% The values %%%
f_0 = figure();
x_axis_in = 0:0.0000625:0.020;
y_input = s(1:321);
hold on
plot(x_axis_in, y_input);
stem(x_axis_dec, y_decoded);
xlabel('t [s]');
ylabel('s[n], symbols');
hold off
saveas(f_0, 'second.png');
close(f_0);
fprintf("Dokončil som dekódovanie symbolov\n");

fprintf("\n3. Úloha\n");
%%% Third task %%%
B = [0.0192 -0.0185 -0.0185 0.0192];
A = [1      -2.8870  2.7997 -0.9113];
%ukazmito(B, A, Fs); %-- kontrola ci je naozaj stabilny
f_1 = figure();
zplane(B, A);
grid
xlabel('Real part');
ylabel('Imaginary part');
saveas(f_1, 'stability.png');
close(f_1);
fprintf("Vytvoril som filter\n");

fprintf("\n4. Úloha\n");
%%% Fourth task %%%
f_2 = figure();
x_axis = (0:255) / 256 * Fs / 2; % -- podla ukazmito
y_axis = abs(freqz(B, A, 256));
plot(x_axis, y_axis);
xlabel('Hz');
saveas(f_2, 'charac.png');
close(f_2);
fprintf("Vytvoril som kmitočtovú charakteristiku filtra\n");

fprintf("\n5. Úloha\n");
%%% Fifth task %%%
f_3 = figure();
ss = filter(B, A, s);
x_axis_s = x_axis_in;
y_axis_s = ss(1:1:321);
hold on
plot(x_axis_in, y_input);
plot(x_axis_s, y_axis_s);
xlabel('t [s]');
ylabel('s[n], ss[n]');
hold off
saveas(f_3, 'filtered.png');
close(f_3);
%% TODO : find out how shifted it is!!!!!
fprintf("Prefiltroval som signál s filtrom\n");

fprintf("\n6. Úloha\n");
%%% Sixth task %%%
f_4 = figure();
s_shifted = ss(1+16:1:32000);
k = 1;
for i=1:period:size(ss)-16
  if(s_shifted(i+8) > 0)
    shift_decoded(k) = 1;
  elseif(s_shifted(i+8) < 0)
    shift_decoded(k) = 0;
  end
  k = k + 1;
end

shift_y_decoded = shift_decoded(1:20);
hold on
plot(x_axis_in, y_input);
plot(x_axis_s, y_axis_s);
stem(x_axis_dec, shift_y_decoded);
plot(x_axis_s, s_shifted(1:1:321));
xlabel('t [s]');
ylabel('s[n], ss[n], ss_{shifted}[n]');
hold off
saveas(f_4, 'shifted.png');
close(f_4);
fprintf("Posunul som signál a všetko zobrazil na grafe\n");

fprintf("\n7. Úloha\n");
%%% Seventh task %%%
incorrect_counter = 0;
size_shifted = length(shift_decoded);

for i=1:1:size_shifted
  if(shift_decoded(i) != decoded(i))
    incorrect_counter = incorrect_counter + 1;
  end
end

fprintf("Počet chýb je %d\n", incorrect_counter);
fprintf("Nesprávnych symbolov je %f percent\n", incorrect_counter/size_shifted*100);

fprintf("\n8. Úloha\n");
%%% Eight task %%%
f_5 = figure();
fft_s = abs(fft(s));
fft_s = abs(fft_s(1:Fs/2+1));
fft_ss = abs(fft(ss));
fft_ss = abs(fft_ss(1:Fs/2+1));
fft_x_axis = (0:Fs/2);

hold on
plot(fft_x_axis, fft_s);
plot(fft_x_axis, fft_ss);
xlabel('[Hz]');
legend('s[n]', 'ss[n]');
hold off
saveas(f_5, 'fourier.png')
close(f_5);
fprintf("Vypočítal som spektrum pomocou DFT\n");

fprintf("\n10. Úloha\n");
%%% Tenth task %%%
f_6 = figure();
k = (-50:50);
R = xcorr(s, 'biased');
R = R(k + length(s));
plot(k, R);
xlabel('k');
saveas(f_6, 'corr.png');
close(f_6);
fprintf("Vypočítal som korelačný koeficient\n");

fprintf("\n11. Úloha\n");
%%% Eleventh task %%%
fprintf('R[0] = %f\n', R(51));
fprintf('R[1] = %f\n', R(52));
fprintf('R[16] = %f\n', R(67));

fprintf("\n12. Úloha\n");
%%% Twelfth task %%%
f_7 = figure(); 
N = length(s);
L = 50;
h = zeros(L, L);
x = linspace(min(s), max(s), 50);
xcol = x(:);
bigx = repmat(xcol, 1, N);
yr = s(:)';
bigy = repmat(yr, L, 1);
[dummy, ind1] = min(abs(bigy - bigx));
ind2 = ind1(2:N);
for ii = 1:N-1,
  d1 = ind1(ii);
  d2 = ind2(ii);
  h(d1, d2) = h(d1, d2) + 1;
end
surf = (x(2) - x(1))^2;
p = h / N / surf;
imagesc(x, x, p);
colorbar;
xlabel('x1');
ylabel('x2');
saveas(f_7, 'distr.png');
close(f_7);
fprintf("Nakreslil som obrázok\n");

fprintf("\n13. Úloha\n");
%%% Thirtienth task %%%
check = sum(sum (p)) * surf; 
disp(['Výsledok integrálu je ' num2str(check)]) 

fprintf("\n14. Úloha\n");
%%% Fourtienth task %%%
printf('R[1] = %f.\n', sum(sum(repmat(x(:), 1, L) .* repmat(x(:)', L, 1) .* p)) * surf);
