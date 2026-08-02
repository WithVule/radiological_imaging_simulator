Q = [0, 0, -100]; % lokacija izvora u cm


% ucitavanje spektara
spek = {
    '60kVp 1mm', '60kVp 17deg 1000Air 0Be 1Al 0Cu 0Sn 0W 0Ta 0Wa.txt';
    '120kVp 1mm', '120kVp 17deg 1000Air 0Be 1Al 0Cu 0Sn 0W 0Ta 0Wa.txt';
    '60kVp 3mm', '60kVp 17deg 1000Air 0Be 3Al 0Cu 0Sn 0W 0Ta 0Wa.txt';
    '120kVp 3mm', '120kVp 17deg 1000Air 0Be 3Al 0Cu 0Sn 0W 0Ta 0Wa.txt'
        };

E = linspace(1, 120, 500)';
N0 = zeros(length(E), 4);

figure;
for i = 1:4
    data_n = readmatrix(spek{i,2}, 'NumHeaderLines', 3);
    En_raw = data_n(:,1);
    N0_raw = data_n(:,2);
    N0(:, i) = interp1(En_raw, N0_raw, E, 'linear', 'extrap');
    N0(:, i) = max(N0(:, i), 0);

    subplot(2, 2, i)
    plot(E, N0(:, i), 'LineWidth', 1.5);
    xlabel('Energija (keV)'); ylabel('N0');
    title(spek(i,1)); grid on;
end


% atenuacija 6. grupe (ime, gustina [g/cm^3], data)
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

% racunanje atenuacija
mu1 = mu(:, 1);
mu2 = mu(:, 2);
mu3 = mu(:, 3);
mu4 = mu(:, 4);
mu5 = mu(:, 5);

%racunanje efikasnosti detektora sa 200um
eps_CsI = 1 - exp(-mu(:, 6) * 200/1e4);
eps_Se = 1 - exp(-mu(:, 7) * 200/1e4);


% detektorski red u cm -- centar u (0, 0, 7)
z_d = 7;
n_d = 380;
l_d = 0.1;

x_d = linspace(-n_d/2*l_d + l_d/2, n_d/2*l_d - l_d/2, n_d);
[X_d, Y_d] = meshgrid(x_d);
centar_piksela = [X_d(:), Y_d(:), z_d*ones(numel(X_d), 1)];

% podrebni racun za duzine
l1 = lenghtTroughElipse(Q, centar_piksela, [0, 0, 0], 30, 30, 6);
l2 = lenghtTroughSphere(Q, centar_piksela, [0, 6, 4], 2);
l3 = lenghtTroughSphere(Q, centar_piksela, [10, 10, 0], 3);
l4 = lenghtTroughSphere(Q, centar_piksela, [-10, -10, 0], 3);
l5 = lenghtTroughCylinder(Q, centar_piksela, [0, 0, -4], 2, 8);


% racun transmitovanog N
N_tr = zeros(size(centar_piksela,1), 4);

for i = 1:length(E)
    N_tr = N_tr + N0(i, :).*exp(-(l1-l2-l3-l4-l5).*mu1(i) -l2.*mu2(i) - l3.*mu3(i) - l4.*mu4(i) - l5.*mu5(i));
end

N_tr_grid = zeros(n_d, n_d, 4);
% prebacivanje iz arraya u matricu
for i = 1:n_d
    for j = 1:n_d
        N_tr_grid(n_d-j+1, i, :) = N_tr((i-1)*n_d + j, :);
    end
end

cmap = gray(256);

figure; title('Idealni detektor');
subplot(2, 2, 1);
render1 = imagesc(N_tr_grid(:, :, 1)); colormap(cmap); colorbar;
title(spek(1, 1));
subplot(2, 2, 2);
render2 = imagesc(N_tr_grid(:, :, 2)); colormap(cmap); colorbar;
title(spek(2, 1));
subplot(2, 2, 3);
render3 = imagesc(N_tr_grid(:, :, 3)); colormap(cmap); colorbar;
title(spek(3, 1));
subplot(2, 2, 4);
render4 = imagesc(N_tr_grid(:, :, 4)); colormap(cmap); colorbar;
title(spek(4, 1));
imcontrast;


% dodatne figure za razlicite efikasnosti detektora

%CsI
N_tr_CsI = zeros(size(centar_piksela,1), 4);

for i = 1:length(E)
    N_tr_CsI = N_tr_CsI + eps_CsI(i) .* N0(i, :) .* exp(-(l1-l2-l3-l4-l5).*mu1(i) -l2.*mu2(i) - l3.*mu3(i) - l4.*mu4(i) - l5.*mu5(i));
end

N_tr_CsI_grid = zeros(n_d, n_d, 4);
% prebacivanje iz arraya u matricu
for i = 1:n_d
    for j = 1:n_d
        N_tr_CsI_grid(n_d-j+1, i, :) = N_tr_CsI((i-1)*n_d + j, :);
    end
end

figure; title('CsI detektor');
subplot(2, 2, 1);
render1_CsI = imagesc(N_tr_CsI_grid(:, :, 1)); colormap(cmap); colorbar;
title(spek(1, 1));
subplot(2, 2, 2);
render2_CsI = imagesc(N_tr_CsI_grid(:, :, 2)); colormap(cmap); colorbar;
title(spek(2, 1));
subplot(2, 2, 3);
render3_CsI = imagesc(N_tr_CsI_grid(:, :, 3)); colormap(cmap); colorbar;
title(spek(3, 1));
subplot(2, 2, 4);
render4_CsI = imagesc(N_tr_CsI_grid(:, :, 4)); colormap(cmap); colorbar;
title(spek(4, 1));
imcontrast;

% Se
N_tr_Se = zeros(size(centar_piksela,1), 4);

for i = 1:length(E)
    N_tr_Se = N_tr_Se + eps_Se(i) .* N0(i, :) .* exp(-(l1-l2-l3-l4-l5).*mu1(i) -l2.*mu2(i) - l3.*mu3(i) - l4.*mu4(i) - l5.*mu5(i));
end

N_tr_Se_grid = zeros(n_d, n_d, 4);
% prebacivanje iz arraya u matricu
for i = 1:n_d
    for j = 1:n_d
        N_tr_Se_grid(n_d-j+1, i, :) = N_tr_Se((i-1)*n_d + j, :);
    end
end

figure; title('CsI detektor');
subplot(2, 2, 1);
render1_Se = imagesc(N_tr_Se_grid(:, :, 1)); colormap(cmap); colorbar;
title(spek(1, 1));
subplot(2, 2, 2);
render2_Se = imagesc(N_tr_Se_grid(:, :, 2)); colormap(cmap); colorbar;
title(spek(2, 1));
subplot(2, 2, 3);
render3_Se = imagesc(N_tr_Se_grid(:, :, 3)); colormap(cmap); colorbar;
title(spek(3, 1));
subplot(2, 2, 4);
render4_Se = imagesc(N_tr_Se_grid(:, :, 4)); colormap(cmap); colorbar;
title(spek(4, 1));
imcontrast;