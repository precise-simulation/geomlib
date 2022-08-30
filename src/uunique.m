function [ b, i, j ] = uunique( a )
%UUNIQUE Unsorted set unique.
%
%   [ B, I, J ] = UUNIQUE( A ) for the vector or row array A returns
%   the same values in the same order as in A without repetitions.
%   I is an index vector so that A = B(I), and J with B = A(J).
%
%   See also UNIQUE

% Initial version 170131.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.


rowvec = 0;
if( any(size(a)==1) )
  a = a(:);
  rowvec = 1;
end
n = size( a, 1 );


[as,is] = sortrows( a );   % Sorting in ascending order, as = a(is).

% Mask indicating duplicated rows.
mask = all( [  repmat(0,1,size(a,2));
              ~[as(2:n,:) - as(1:n-1,:)] ], 2 );


% Compute start 's' and end 'e' position of duplicate rows,
% and the number of duplicates for each group 'ndup'.
g  = diff( [0; mask; 0] );   % Grouping, 1 indicates start of group, -1 position after end of group.
s  = find( g== 1 );
e  = find( g==-1 ) - 1;
r  = is( s - 1 );   % Original position of reference value for each group (duplicate to keep).
ndup = e - s + 1;

% Index positions to duplicates originals.
ix([cumsum([1; ndup(ndup>0)])]) = 1;
idup = r( cumsum( ix(1:find(ix,1,'last')-1) ) );

% Original position of duplicates.
% isdi = cell2mat( arrayfun( @colon, s', e', 'UniformOutput', false ) )';
isdi = zeros( max(e)+1, 1 );
isdi(s) = 1;
isdi(e+1) = isdi(e+1) - 1;
isdi = find( cumsum(isdi) );

% Compute index order without duplicates.
i = zeros( n, 1 );
j = sort(is(mask==0));
i(j) = 1:sum(mask==0);

% Assign index to duplicates.
i(is(isdi)) = i(idup);


% Output.
if( rowvec )
  i = i';
  b = a(j)';
else
  b = a(j,:);
end



%% UNIT TESTS

%!assert( isequal( uunique(1:10),    1:10 ) )
%!assert( isequal( uunique([1:10]'), 1:10 ) )
%!assert( isequal( uunique(10:-2:1), 10:-2:1 ) )

%!test % Test 1
%! a = [ 11 18 12 14 13 15 13 16 17 18 18 ];
%! [b,i] = uunique(a);
%! assert( isequal( b, [11 18 12 14 13 15 16 17] ) )
%! assert( isequal( i, [1 2 3 4 5 6 5 7 8 2 2] ) )
%! assert( isequal( a, b(i) ) )
%! [b,i] = uunique(a');
%! assert( isequal( b, [11 18 12 14 13 15 16 17] ) )
%! assert( isequal( i, [1 2 3 4 5 6 5 7 8 2 2] ) )
%! assert( isequal( a, b(i) ) )

%!test % Test 2
%! a = repmat( [ 11 18 12 14 13 15 13 16 17 18 18 ]', 1, 2 );
%! [b,i] = uunique( a );
%! assert( isequal( b, repmat([11 18 12 14 13 15 16 17]',1,2) ) )
%! assert( isequal( i, [1 2 3 4 5 6 5 7 8 2 2]' ) )
%! assert( isequal( a, b(i,:) ) )

%!test % Test 3
%! a = repmat( [ 11 18 12 14 13 15 13 16 17 18 18 ]', 1, 3 );
%! [b,i] = uunique( a );
%! assert( isequal( b, repmat([11 18 12 14 13 15 16 17]',1,3) ) )
%! assert( isequal( i, [1 2 3 4 5 6 5 7 8 2 2]' ) )
%! assert( isequal( a, b(i,:) ) )

%!test % Test 4
%! a = [ 13 13 6 17 17 17 ];
%! [b,i] = uunique(a);
%! assert( isequal( b, [13 6 17]) )
%! assert( isequal( i, [1 1 2 3 3 3] ) )
%! assert( isequal( a, b(i) ) )
