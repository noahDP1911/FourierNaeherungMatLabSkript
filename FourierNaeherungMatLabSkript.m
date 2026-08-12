% Parameter
A = 1;              % Amplitude für Bit 1
A_0 = 0;           % Amplitude für Bit 0
T_bit = 1;          % Dauer eines Bits
bit_pattern = [1 1 0 1 0 0 0 1]; % Zu übertragendes Bitmuster
T_total = T_bit * length(bit_pattern); % Gesamtdauer einer Periode
f1 = 1 / T_total; % Grundfrequenz des Signals

% Zeitvektor mit 100 Abtastpunkten pro Bit und dem Periodenende
t = linspace(0, T_total, 100 * length(bit_pattern) + 1);

% Rechtecksignal erzeugen basierend auf Bitmuster
x_t = zeros(size(t));
for i = 1:length(bit_pattern)
    % Start- und Endzeit des aktuellen Bits
    t_start = (i - 1) * T_bit;
    t_end = i * T_bit;

    % Setze Amplitude für das aktuelle Bit
    if bit_pattern(i) == 1
        x_t(t >= t_start & t < t_end) = A;
    else
        x_t(t >= t_start & t < t_end) = A_0;
    end
end

% Der letzte Zeitpunkt entspricht dem Beginn der nächsten Periode.
x_t(end) = x_t(1);

% Anzahl der Harmonischen N
N = 5; % Anzahl der zu berücksichtigenden Harmonischen

% Gleichanteil berechnen und Speicher für die Harmonischen anlegen
a_0 = (2 / T_total) * trapz(t, x_t);
x_approx = (a_0 / 2) * ones(size(t));
harmonics = zeros(N, length(t));
for n = 1:N
    % Berechne den Fourier-Koeffizienten (a_n und b_n)
    omega_n = 2 * pi * n * f1;
    a_n = (2 / T_total) * trapz(t, x_t .* cos(omega_n * t));
    b_n = (2 / T_total) * trapz(t, x_t .* sin(omega_n * t));

    % Speichere und addiere den Beitrag der n-ten Harmonischen
    harmonics(n, :) = a_n * cos(omega_n * t) + b_n * sin(omega_n * t);
    x_approx = x_approx + harmonics(n, :);
end

% Berechne den Fehler (RMS-Fehler)
rms_error = sqrt(mean((x_t - x_approx).^2));

% Anzeige des Fehlers
disp(['RMS-Fehler mit ', num2str(N), ' Harmonischen: ', num2str(rms_error)])

% Plot: Originales Rechtecksignal und Näherung
figure;
subplot(3, 1, 1);
plot(t, x_t, 'LineWidth', 2);
title(['Originales Rechtecksignal (Bitmuster: ', num2str(bit_pattern), ')']);
xlabel('Zeit (s)');
ylabel('Amplitude');
grid on;

subplot(3, 1, 2);
plot(t, harmonics.', 'LineWidth', 1);
title(['Einzelbeiträge der ersten ', num2str(N), ' Harmonischen']);
xlabel('Zeit (s)');
ylabel('Amplitude');
grid on;

subplot(3, 1, 3);
plot(t, x_approx, 'LineWidth', 2);
title(['Fourier-Näherung mit ', num2str(N), ' Harmonischen']);
xlabel('Zeit (s)');
ylabel('Amplitude');
grid on;
