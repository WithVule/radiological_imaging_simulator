mat = {'Gd2O2S', 7.32, 'Gd2O2S.txt';
    'BaFBr',  4.90, 'BaFBr.txt';
    'CsI',    4.51, 'CsI.txt';
    'Se',     4.81, 'Se.txt'};

E = linspace(1, 140, 500)';
mu = zeros(length(E), 4);

for i = 1:4
    data = readmatrix(mat{i,3}, 'NumHeaderLines', 3);
    E_raw = data(:,1) * 1000;
    mu_rho_raw = data(:,2);

    % K ivica duplikati fix
    for j = 2:length(E_raw)
        if E_raw(j) == E_raw(j-1)
            E_raw(j) = E_raw(j-1) + eps(E_raw(j-1));
        end
    end

    mu_rho = interp1(E_raw, mu_rho_raw, E, 'linear', 'extrap');
    mu(:,i) = mu_rho * mat{i,2};
end

eps0 = 1 - exp(-mu * 200/1e4);
figure; plot(E, eps0, 'LineWidth', 1.5);
xlabel('Energija (keV)'); ylabel('\epsilon_i');
title('200 μm debljina'); legend(mat(:,1)); grid on;

d = [200; zeros(3,1)];
eps_ref = 1 - exp(-mu(:,1) * 200/1e4);

for i = 2:4
    d_cm = median(-log(1 - eps_ref) ./ mu(:,i), 'omitnan');
    d(i) = d_cm * 1e4;
    eps_ref = 1 - exp(-mu(:,i) * d_cm);
end

eps_opt = 1 - exp(-mu * (d ./ 1e4));

fprintf('Debljine(μm):\n');
fprintf('%s: %.0f\n', mat{1,1}, d(1));
fprintf('%s: %.0f\n', mat{2,1}, d(2));
fprintf('%s: %.0f\n', mat{3,1}, d(3));
fprintf('%s: %.0f\n', mat{4,1}, d(4));
