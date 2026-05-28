Q = [0, 0, -100]; % lokacija izvora u cm
Eq = [60, 120]; % energije izvora u keV

% anoda od W

d_filt = [0.1, 0.3]; % debljina filtra u cm
an_ugao = 17; % anodni ugao

% dimenzije elipsoida - fantoma u cm -- centar u (0, 0, 0)
a = 30; %x
b = 30; %y
c = 6;  %z

function [z1, z2] = elipsoid(x, y)
    z1 = c*sqrt(1-(x/a)^2 - (y/b)^2);
    z2 = -z1;
end

