function [scattered, D] = scatterCheck(xi, Q, xyz, mu, E_grid, E_ph, rho)

% podrebni racun za duzine
l1 = lenghtTroughElipse(Q, xyz, [0, 0, 0], 30, 30, 6);
l2 = lenghtTroughSphere(Q, xyz, [0, 6, 4], 2);
l3 = lenghtTroughSphere(Q, xyz, [10, 10, 0], 3);
l4 = lenghtTroughSphere(Q, xyz, [-10, -10, 0], 3);
l5 = lenghtTroughCylinder(Q, xyz, [0, 0, -4], 2, 8);

mu_E = interp1(E_grid, mu(:,1:5), E_ph, 'linear', 'extrap');

d_scatter = - log(xi);

z1 = mu_E(:,1).*l1;
z2 = z1 + mu_E(:,2).*l2;
z3 = z2 + mu_E(:,3).*l3;
z4 = z3 + mu_E(:,4).*l4;
z5 = z3 + mu_E(:,5).*l5;

scattered = z5 > d_scatter;

sloj = zeros(size(E_ph));
sloj(scattered & d_scatter <= z1) = 1;
sloj(scattered & d_scatter > z1  & d_scatter <= z2) = 2;
sloj(scattered & d_scatter > z2  & d_scatter <= z3) = 3;
sloj(scattered & d_scatter > z3  & d_scatter <= z4) = 4;
sloj(scattered & d_scatter > z4  & d_scatter <= z5) = 5;

Di = zeros(size(E_ph));
for m = 1:5
    idx = sloj == m;
    Di(idx) = E_ph(idx) / rho(m);
end

D = Di(2:5);

end