h = [0.1, 0.05, 0.01, 0.005, 0.001]';
x0 = ones(5,1);
u = sin(x0);

du_plus = (sin(x0 + h) - sin(x0))./h;
du_moins = (sin(x0) - sin(x0-h))./h;
du_cent = 1/2 * (du_plus + du_moins);
du_v = cos(x0);

format shorte

[h, du_v - du_plus, du_v - du_moins, du_v - du_cent]

format longg