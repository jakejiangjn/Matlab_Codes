function [ A_index B_index ] = matchsqn( vectorA, vectorB, radius )
% Find the common or nearest scalar elements in 2 vectors
% Only if these 2 vectors are pre-sorted
% Return the index
if length( vectorA ) < length( vectorB )
    a = vectorA;    b = vectorB;    flag = 1;
else
    a = vectorB;    b = vectorA;    flag = 0;
end;
N = length(a);
A_index = zeros( N, 1 );    B_index = zeros( N, 1 );
index = 1;
for cnt = 1:N
    tmp = abs( b - a(cnt) );
    index_tmp = find( tmp == min(tmp), 1 );
    if tmp(index_tmp) > radius
        continue;
    end;
    A_index(index) = cnt;
    B_index(index) = index_tmp;
    index = index + 1;
end;
A_index = A_index( A_index > 0 );
B_index = B_index( B_index > 0 );
if flag == 0
    tmp = A_index;
    A_index = B_index;
    B_index = tmp;
end;