clc;

tic

max_magic = 1e5;
% Collect the results as they become available.
magicResults = cell(1,max_magic);
for idx = 1:max_magic
  magicResults{idx} = magic(10);  
end
fprintf('Got result...');

toc