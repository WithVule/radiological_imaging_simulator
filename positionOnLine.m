function [x_photon, y_photon] = positionOnLine(xi, Q, z_d)

theta = acos(1-2*xi(:,1));
phi = 2*pi*xi(:,2);

t = (z_d - Q(:,3))./cos(theta);

x_photon = Q(:,1) + sin(theta).*cos(phi).*t;
y_photon = Q(:,2) + sin(theta).*sin(phi).*t;

end