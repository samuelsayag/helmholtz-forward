clc; clear variables; close all;

tic


p = gcp();
max_magic = 1e5 ;
step = 1e4;
% step = fix(max_magic/10);
start = step;
myrem = mod(max_magic, step);

if (myrem)
    v_ind = [start:step:max_magic, max_magic];
else
    v_ind = start:step:max_magic;    
end
    



N = size(v_ind, 2);
% To request multiple evaluations, use a loop.
func = @(x,y) create_magic(x,y);
begi = 1;
for idx = 1:N    
    endi = v_ind(idx);
    f(idx) = parfeval(p, func, 3, begi, endi); % Square size determined by idx
    begi = begi + step;
end
% Collect the results as they become available.
magicResults = cell(1,max_magic);
for idx = 1:N
    % fetchNext blocks until next results are available.
    [completedIdx, value, begi, endi] = fetchNext(f);
    magicResults(begi:endi) = value;
%     fprintf('Got result with index: %d.\n', completedIdx);
end

toc

