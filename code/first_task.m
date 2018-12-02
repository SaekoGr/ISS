%%% First task %%%
[s, Fs] = audioread('xgregu02.wav');
period = 16;  % this is one period
%fileID = fopen('tmp_xgregu02.txt', 'w');
k = 1;
i = 9;

fprintf("Vzorkovacia frekvencia je: %f [Hz].\n", Fs);
fprintf("Dĺžka signálu vo vzorkoch: %f, v sekundách: %f [s].\n", length(s), length(s)/Fs);

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
x_axis_in = 0:0.0000625:0.020;
y_input = s(1:321);
f_0 = figure();
hold on
plot(x_axis_in, y_input);
stem(x_axis_dec, y_decoded);
xlabel('t [s]');
ylabel('s[n], symbols');
hold off
saveas(f_0, 'second.png');
close(f_0);

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

%%% Fourth task %%%
f_2 = figure();
x_axis = (0:255) / 256 * Fs / 2; % -- podla ukazmito
y_axis = abs(freqz(B, A, 256));
plot(x_axis, y_axis);
xlabel('Hz');
saveas(f_2, 'charac.png');
close(f_2);

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

%%% Sixth task %%%
f_4 = figure();
s_shifted = ss(1+8:1:321+8);
hold on
plot(x_axis_in, y_input);
plot(x_axis_s, y_axis_s);
plot(x_axis_s, s_shifted);
xlabel('t [s]');
ylabel('s[n], ss[n], ss_{shifted}[n]');
hold off
saveas(f_4, 'shifted.png');
close(f_4);

%%% Tenth task %%%
f_6 = figure();
k = (-50:50);
R = xcorr(s, 'biased');
R = R(k + length(s));
plot(k, R);
xlabel('k');
saveas(f_6, 'corr.png');
close(f_6);

%%% Eleventh task %%%
fprintf('R[0] = %f\n', R(51));
fprintf('R[1] = %f\n', R(52));
fprintf('R[16] = %f\n', R(67));

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

%%% Thirtienth task %%%
check = sum(sum (p)) * surf; 
disp(['hist2: check -- 2d integral should be 1 and is ' num2str(check)]) 

%%% Fourtienth task %%%
printf('R[1] = %f.\n', sum(sum(repmat(x(:), 1, L) .* repmat(x(:)', L, 1) .* p)) * surf);
