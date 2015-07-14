function [ res, begi, endi ] = create_magic( begi, endi )
%CREATE_MAGIC Summary of this function goes here
%   Detailed explanation goes here



res = cell(1, endi - begi + 1);

for idx = 1: size(res,2)
  res{idx} = magic(10);  
end

end

