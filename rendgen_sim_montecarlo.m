clear; clc; close all;

% ODREDJIVANJE POZICIJE IZVORA
txt_source = input("Da li koristim model žižne tačke? Y/N [Y]: ", "s");
if isempty(txt_source)
    txt_source = 'Y';
end

if txt_source == 'N'
    point_source = 1;
elseif txt_source == 'Y'
    point_source = 0;
else
    printf("Nevažeći odgovor.");
    quit;
end
% U slucaju da je izabran model zizne tacke, onda se simulira za svaki
% foton izlazna tacka, woohoo

% broj fotona
N_ph = 10^8;

% definisanje potrebnih matrica
photons = zeros(N_ph, 1);
x_Q = zeros(N_ph, 1);
y_Q = zeros(N_ph, 1);
x_d = zeros(N_ph, 1);
y_d = zeros(N_ph, 1);
d_ph = zeros(N_ph, 1);

% x_Q, y_Q, energy, theta, phi, attenuation
R = rand(N_ph, 6);

Q = zeros(N_ph, 3);

if ~point_source
    dim = [0.1, 0.2]; % u cm
    [x_Q, y_Q] = positionOnFocalPoint(R(:, 1:2), dim);
    Q(:, 1) = x_Q;
    Q(:, 2) = y_Q;
    Q(:, 3) = -100; % pozicija izvora u cm
else
    Q = repmat([0, 0, -100], N_ph, 1);
end


% ODREDJIVANJE ENERGIJE FOTONA
spectras = {
    '60kVp 1mm', '60kVp 17deg 1000Air 0Be 1Al 0Cu 0Sn 0W 0Ta 0Wa.txt';
    '120kVp 1mm', '120kVp 17deg 1000Air 0Be 1Al 0Cu 0Sn 0W 0Ta 0Wa.txt';
    '60kVp 3mm', '60kVp 17deg 1000Air 0Be 3Al 0Cu 0Sn 0W 0Ta 0Wa.txt';
    '120kVp 3mm', '120kVp 17deg 1000Air 0Be 3Al 0Cu 0Sn 0W 0Ta 0Wa.txt'
    };

E = linspace(1, 120, 500)';
N0 = zeros(length(E), 4);

figure;
for i = 1:4
    data_n = readmatrix(spectras{i,2}, 'NumHeaderLines', 3);
    En_raw = data_n(:,1);
    N0_raw = data_n(:,2);
    N0(:, i) = interp1(En_raw, N0_raw, E, 'linear', 'extrap');
    N0(:, i) = max(N0(:, i), 0);

    subplot(2, 2, i)
    plot(E, N0(:, i), 'LineWidth', 1.5);
    xlabel('Energija (keV)'); ylabel('N0');
    title(string(i) + " - "+ spectras(i,1)); grid on;
end

n_spectra = input("Koji spektar energije koristim? 1-4 [1]: ");
if isempty(n_spectra)
    n_spectra = 1;
end

E_ph = photonEnergy(R(:, 3), E, N0(:, n_spectra));


% ODREĐIVANJE LIN. KOEF. ATTENUACIJE
materijali = {
    'H2O', 1, 'H2O.txt';
    'Ovary', 1.050, 'Ovary.txt';
    'Bone, Cortical', 1.920, 'Bone.txt';
    'Brain, Gray/White Matter', 1.040, 'Brain.txt';
    'Adipose Tissue', 0.950, 'Adipose.txt'
    'CsI', 4.51, 'CsI.txt';
    'Se', 4.81, 'Se.txt'
    };

mu = zeros(length(E), 7);

for i = 1:7
    data_mu = readmatrix(materijali{i,3}, 'NumHeaderLines', 3);
    Emu_raw = data_mu(:,1) * 1000;
    mu_rho_raw = data_mu(:,2);

    % K ivica duplikati fix
    for j = 2:length(Emu_raw)
        if Emu_raw(j) == Emu_raw(j-1)
            Emu_raw(j) = Emu_raw(j-1) + eps(Emu_raw(j-1));
        end
    end

    mu_rho = interp1(Emu_raw, mu_rho_raw, E, 'linear', 'extrap');
    mu(:,i) = mu_rho * materijali{i,2};
end

%racunanje efikasnosti detektora sa 200um
eps_CsI = 1 - exp(-mu(:, 6) * 200/1e4);
eps_Se = 1 - exp(-mu(:, 7) * 200/1e4);

% detektorski red u cm -- centar u (0, 0, 7)
z_d = 7; % promenljivo
n_d = 380;
l_d = 0.1;

[x_photon, y_photon] = positionOnLine(R(:,4:5), Q, z_d);
z_photon = z_d*ones(size(x_photon));

[ix, iy, valid] = pixelIndex(x_photon, y_photon, n_d, l_d);

[scattered, D] = scatterCheck(R(valid,6), Q(valid,:), [x_photon(valid), y_photon(valid), z_photon(valid)], mu(:,1:5), E, E_ph(valid), [materijali{1:5,2}]);

not_scattered_full = ones(N_ph, 1);
not_scattered_full(valid) = not_scattered_full(valid) - scattered;

final_mask = valid & not_scattered_full;

E_ph_final = E_ph(valid);

slika_ideal = weightedImage(ix, iy, final_mask, E_ph, E, 1, n_d);
slika_CsI = weightedImage(ix, iy, final_mask, E_ph, E, eps_CsI, n_d);
slika_Se = weightedImage(ix, iy, final_mask, E_ph, E, eps_Se, n_d);

cmap = gray(256);
figure;
subplot(1,3,1); imagesc(slika_ideal); colormap(cmap); colorbar; axis image;
title('Idealan detektor');

subplot(1,3,2); imagesc(slika_CsI); colormap(cmap); colorbar; axis image;
title('CsI, 200mm');

subplot(1,3,3); imagesc(slika_Se); colormap(cmap); colorbar; axis image;
title('Se, 200mm');

disp(D);