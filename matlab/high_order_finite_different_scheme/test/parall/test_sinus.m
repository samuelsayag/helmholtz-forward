wclc;
tic

sof = 2.^16

% for i = 1:sof
parfor i = 1:sof
  A(i) = sin(i*2*pi/sof);
end
plot(A)

toc