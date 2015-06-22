mygrid = linspace(0.125,0.375,3);
% mygrid = linspace(0,0.5,3);
x = repmat(mygrid,3,1);
y = repmat(flip(mygrid'),1,3);

create_grid_figure(x,y)


