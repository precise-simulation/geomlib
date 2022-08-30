function [ res ] = test_csg_op0( fun )
if(nargin),if(ischar(fun)),fun=str2func(fun);end,fun();if(nargout>0),res=[];end,return,end
res = mloct_test_caller(cellfun(@str2func,get_subfuns([mfilename('fullpath'),'.m'],'test'),'UniformOutput',0));


function test_build_extract_usquare
% Test 1 - Build and extract unit square.

rectangle = { 1, 2, 3, 4 ;
              2, 2, 2, 2 ;
              [0,-1], [1,0], [0,1], [-1,0] ;
              [0,0;1,0], [1,0;1,1], [1,1;0,1], [0,1;0,0] ;
              [1/2,0], [1,1/2], [1/2,1], [0,1/2] ;
              1, 1, 1, 1 };

a = csg_op( rectangle, 'b' );

a_ref = { [ 0, 0 ], [ 0, -1 ], [], { [ 1, 0 ], [ 1, 0 ], [], { [ 1, 1], ...
          [ 0, 1 ], [], { [ 0, 1 ], [ -1, 0 ], [], [], { 4; 2; [ -1, 0 ];
          [ 0, 1; 0, 0 ]; [ 0, 0.5 ]; 1 } }, { 3; 2; [ 0, 1 ]; [ 1, 1; 0, 1 ];
          [ 0.5, 1 ]; 1 } }, { 2; 2; [ 1, 0 ]; [ 1, 0; 1, 1 ]; [ 1, 0.5 ]; 1 } }, ...
          { 1; 2; [ 0, -1 ]; [ 0, 0; 1, 0 ]; [ 0.5, 0 ]; 1 } };

assert( isequal(a,a_ref) )

polygons = csg_op( a, 'e' );

assert( isequal(rectangle,polygons) )


function test_build_extract_ucube
% Test 2 - Build and extract unit cube.

block = { 1, 2, 3, 4, 5, 6 ;
          4, 4, 4, 4, 4, 4 ;
          [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
          [0,0,0;0,1,0;1,1,0;1,0,0], [0,0,0;1,0,0;1,0,1;0,0,1], [1,0,0;1,1,0;1,1,1;1,0,1], ...
          [1,1,0;0,1,0;0,1,1;1,1,1], [0,1,0;0,0,0;0,0,1;0,1,1], [0,0,1;1,0,1;1,1,1;0,1,1] ;
          [1/2,1/2,0], [1/2,0,1/2], [1,1/2,1/2], [1/2,1,1/2], [0,1/2,1/2], [0,1/2,1/2] ;
          1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

a = csg_op( block, 'b' );

a_ref = { [ 0, 0, 0 ], [ 0, 0, -1 ], [] , { [ 0, 0, 0 ], [ 0, -1, 0 ], ...
          [] , { [ 1, 0, 0 ], [ 1, 0, 0 ], [] , { [ 1, 1, 0 ], [ 0, 1, 0 ], [] , ...
          { [ 0, 1, 0 ], [ -1, 0, 0 ], [], { [ 0, 0, 1 ], [ 0, 0, 1 ], [], [], ...
          { 6; 4; [ 0, 0, 1 ]; [ 0, 0, 1; 1, 0, 1; 1, 1, 1; 0, 1, 1 ]; [ 0, 0.5, ...
          0.5 ]; 1/sqrt(2) } }, { 5; 4; [ -1, 0, 0 ]; [ 0, 1, 0; 0, 0, 0; 0, 0, 1;
          0, 1, 1 ]; [ 0, 0.5, 0.5 ]; 1/sqrt(2) } }, { 4; 4; [ 0, 1, 0 ]; [ 1, 1, ...
          0; 0, 1, 0; 0, 1, 1; 1, 1, 1 ]; [ 0.5, 1, 0.5 ]; 1/sqrt(2) } }, { 3; 4;
          [ 1, 0, 0 ]; [ 1, 0, 0; 1, 1, 0; 1, 1, 1; 1, 0, 1 ]; [ 1, 0.5, 0.5 ];
          1/sqrt(2) } }, { 2; 4; [ 0, -1, 0 ]; [ 0, 0, 0; 1, 0, 0; 1, 0, 1; 0, 0, ...
          1 ]; [ 0.5, 0, 0.5 ]; 1/sqrt(2) } }, { 1; 4; [ 0, 0, -1 ]; [ 0, 0, 0; 0, ...
          1, 0; 1, 1, 0; 1, 0, 0 ]; [ 0.5, 0.5, 0 ]; 1/sqrt(2) } };

assert( isequal(a,a_ref) )

polygons = csg_op( a, 'e' );

assert( isequal(block,polygons) )


function test_join_usquares
% Test 3 - Join two unit squares.

r1 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0], [1,0;1,1], [1,1;0,1], [0,1;0,0] ;
       [1/2,0], [1,1/2], [1/2,1], [0,1/2] ;
       1, 1, 1, 1 };

a = csg_op( r1, 'b' );

r2 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0]+1/2, [1,0;1,1]+1/2, [1,1;0,1]+1/2, [0,1;0,0]+1/2 ;
       [1/2,0]+1/2, [1,1/2]+1/2, [1/2,1]+1/2, [0,1/2]+1/2 ;
       1, 1, 1, 1 };

b = csg_op( r2, 'b' );

c = csg_op( a, b, '+' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)>=8 && size(polygons,2)<=10 )

[tmp1,tmp2,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( all(n_cnt==2) )


function test_subtract_usquares1
% Test 4 - Subtract r1 - r2.

r1 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0], [1,0;1,1], [1,1;0,1], [0,1;0,0] ;
       [1/2,0], [1,1/2], [1/2,1], [0,1/2] ;
       1, 1, 1, 1 };

a = csg_op( r1, 'b' );

r2 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0]+1/2, [1,0;1,1]+1/2, [1,1;0,1]+1/2, [0,1;0,0]+1/2 ;
       [1/2,0]+1/2, [1,1/2]+1/2, [1/2,1]+1/2, [0,1/2]+1/2 ;
       1, 1, 1, 1 };

b = csg_op( r2, 'b' );

c = csg_op( a, b, '-' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)==6 )

[tmp1,tmp2,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( all(n_cnt==2) )


function test_subtract_usquares2
% Test 5 - Subtract r2 - r1.

r1 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0], [1,0;1,1], [1,1;0,1], [0,1;0,0] ;
       [1/2,0], [1,1/2], [1/2,1], [0,1/2] ;
       1, 1, 1, 1 };

a = csg_op( r1, 'b' );

r2 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0]+1/2, [1,0;1,1]+1/2, [1,1;0,1]+1/2, [0,1;0,0]+1/2 ;
       [1/2,0]+1/2, [1,1/2]+1/2, [1/2,1]+1/2, [0,1/2]+1/2 ;
       1, 1, 1, 1 };

b = csg_op( r2, 'b' );

c = csg_op( b, a, '-' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)==6 )

[tmp1,tmp2,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( all(n_cnt==2) )


function test_intersect_usquares
% Test 6 - Inersect r2 & r1.

r1 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0], [1,0;1,1], [1,1;0,1], [0,1;0,0] ;
       [1/2,0], [1,1/2], [1/2,1], [0,1/2] ;
       1, 1, 1, 1 };

a = csg_op( r1, 'b' );

r2 = { 1, 2, 3, 4 ;
       2, 2, 2, 2 ;
       [0,-1], [1,0], [0,1], [-1,0] ;
       [0,0;1,0]+1/2, [1,0;1,1]+1/2, [1,1;0,1]+1/2, [0,1;0,0]+1/2 ;
       [1/2,0]+1/2, [1,1/2]+1/2, [1/2,1]+1/2, [0,1/2]+1/2 ;
       1, 1, 1, 1 };

b = csg_op( r2, 'b' );

c = csg_op( b, a, '&' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)==4 )

[tmp1,tmp2,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( all(n_cnt==2) )


function test_join_ucubes
% Test 7 - Join two unit cubes.
b1 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0], [0,0,0;1,0,0;1,0,1;0,0,1], [1,0,0;1,1,0;1,1,1;1,0,1], ...
       [1,1,0;0,1,0;0,1,1;1,1,1], [0,1,0;0,0,0;0,0,1;0,1,1], [0,0,1;1,0,1;1,1,1;0,1,1] ;
       [1/2,1/2,0], [1/2,0,1/2], [1,1/2,1/2], [1/2,1,1/2], [0,1/2,1/2], [0,1/2,1/2] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

a = csg_op( b1, 'b' );

b2 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0]+1/2, [0,0,0;1,0,0;1,0,1;0,0,1]+1/2, [1,0,0;1,1,0;1,1,1;1,0,1]+1/2, ...
       [1,1,0;0,1,0;0,1,1;1,1,1]+1/2, [0,1,0;0,0,0;0,0,1;0,1,1]+1/2, [0,0,1;1,0,1;1,1,1;0,1,1]+1/2 ;
       [1,1,1/2], [1,1/2,1], [3/2,1,1], [1,3/2,1], [1/2,1,1], [1/2,1,1] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

b = csg_op( b2, 'b' );

c = csg_op( a, b, '+' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)>=9 && size(polygons,2)<=21 )

[tmp1,iind,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( length(iind)==28 )
assert( all(n_cnt>=2 & n_cnt<=4) )


function test_subtract_ucubes1
% Test 8 - Subtract b1 - b2.

b1 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0], [0,0,0;1,0,0;1,0,1;0,0,1], [1,0,0;1,1,0;1,1,1;1,0,1], ...
       [1,1,0;0,1,0;0,1,1;1,1,1], [0,1,0;0,0,0;0,0,1;0,1,1], [0,0,1;1,0,1;1,1,1;0,1,1] ;
       [1/2,1/2,0], [1/2,0,1/2], [1,1/2,1/2], [1/2,1,1/2], [0,1/2,1/2], [0,1/2,1/2] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

a = csg_op( b1, 'b' );

b2 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0]+1/2, [0,0,0;1,0,0;1,0,1;0,0,1]+1/2, [1,0,0;1,1,0;1,1,1;1,0,1]+1/2, ...
       [1,1,0;0,1,0;0,1,1;1,1,1]+1/2, [0,1,0;0,0,0;0,0,1;0,1,1]+1/2, [0,0,1;1,0,1;1,1,1;0,1,1]+1/2 ;
       [1,1,1/2], [1,1/2,1], [3/2,1,1], [1,3/2,1], [1/2,1,1], [1/2,1,1] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

b = csg_op( b2, 'b' );

c = csg_op( a, b, '-' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)==12 )

[tmp1,tmp2,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( all(n_cnt>=2 & n_cnt<=3) )


function test_subtract_ucubes2
% Test 9 - Subtract b2 - b1.

b1 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0], [0,0,0;1,0,0;1,0,1;0,0,1], [1,0,0;1,1,0;1,1,1;1,0,1], ...
       [1,1,0;0,1,0;0,1,1;1,1,1], [0,1,0;0,0,0;0,0,1;0,1,1], [0,0,1;1,0,1;1,1,1;0,1,1] ;
       [1/2,1/2,0], [1/2,0,1/2], [1,1/2,1/2], [1/2,1,1/2], [0,1/2,1/2], [0,1/2,1/2] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

a = csg_op( b1, 'b' );

b2 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0]+1/2, [0,0,0;1,0,0;1,0,1;0,0,1]+1/2, [1,0,0;1,1,0;1,1,1;1,0,1]+1/2, ...
       [1,1,0;0,1,0;0,1,1;1,1,1]+1/2, [0,1,0;0,0,0;0,0,1;0,1,1]+1/2, [0,0,1;1,0,1;1,1,1;0,1,1]+1/2 ;
       [1,1,1/2], [1,1/2,1], [3/2,1,1], [1,3/2,1], [1/2,1,1], [1/2,1,1] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

b = csg_op( b2, 'b' );

c = csg_op( b, a, '-' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)==12 )

[tmp1,tmp2,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( all(n_cnt>=2 & n_cnt<=3) )


function test_intersect_ucubes
% Test 10 - Intersect b1 & b2.

b1 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0], [0,0,0;1,0,0;1,0,1;0,0,1], [1,0,0;1,1,0;1,1,1;1,0,1], ...
       [1,1,0;0,1,0;0,1,1;1,1,1], [0,1,0;0,0,0;0,0,1;0,1,1], [0,0,1;1,0,1;1,1,1;0,1,1] ;
       [1/2,1/2,0], [1/2,0,1/2], [1,1/2,1/2], [1/2,1,1/2], [0,1/2,1/2], [0,1/2,1/2] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

a = csg_op( b1, 'b' );

b2 = { 1, 2, 3, 4, 5, 6 ;
       4, 4, 4, 4, 4, 4 ;
       [0,0,-1], [0,-1,0], [1,0,0], [0,1,0], [-1,0,0], [0,0,1] ;
       [0,0,0;0,1,0;1,1,0;1,0,0]+1/2, [0,0,0;1,0,0;1,0,1;0,0,1]+1/2, [1,0,0;1,1,0;1,1,1;1,0,1]+1/2, ...
       [1,1,0;0,1,0;0,1,1;1,1,1]+1/2, [0,1,0;0,0,0;0,0,1;0,1,1]+1/2, [0,0,1;1,0,1;1,1,1;0,1,1]+1/2 ;
       [1,1,1/2], [1,1/2,1], [3/2,1,1], [1,3/2,1], [1/2,1,1], [1/2,1,1] ;
       1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2) };

b = csg_op( b2, 'b' );

c = csg_op( a, b, '&' );

polygons = csg_op( c, 'e' );

assert( size(polygons,2)==6 )

[tmp1,tmp2,jind] = l_deduplicate( vertcat(polygons{4,:}), 1 );
n_cnt = hist( jind, 1:max(jind) );

assert( isequal(n_cnt,3*ones(1,8)) )


function [ b, i, j ] = l_deduplicate( a, dim, tol, mode )
%DEDUPLICATE Remove duplicate rows or columns within tolerance.
%
%   [ B, I, J ] = DEDUPLICATE( A, DIM, TOL, MODE ) Removes duplicate rows
%   or columns (according to DIM) from array A with tolerance TOL (default
%   eps*1e3). Returns an array B and optionally index vectors I and J so
%   that B = A(I,:) and A = B(J,:). (Retains the order of the unique
%   rows/columns in A). Additionally MODE is a flag to indicate the use
%   of unique built-in function (MODE = 0) or not (default).
%
%   See also UNIQUE

% Initial version 160513.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.


USE_UNIQUE = 0;

if( nargin>0 )
  if( nargin<4 )
    mode = 1;
    if( nargin<3 )
      tol = eps*1e3;
      if( nargin<2 )
        dim = 1;
      end,end,end
else
  if( nargout==0 )
    help deduplicate, return
  else
    error( 'DEDUPLICATE: Not enough input arguments.' );
  end
end
if( isempty(a) )
  b = [];
  i = [];
  j = [];
  return
end


if( dim==2 )
  a = a';   % Transpose in order to work row wise.
end

if( tol>0 )
  span = tol*max( max(a) - min(a) );
else
  span = abs( tol );
end
spaninv = 1/span;


if( mode==0 )

  [~,i,j] = unique( span*round( spaninv*a ), 'rows', 'first' );

else

  [c,k] = sortrows( span*round( spaninv*a ) );
  n_rows = size( c, 1 );
  not_match = any( c(1:n_rows-1,:)~=c(2:n_rows,:), 2 );
  % j = k;
  j(k) = cumsum( [ 1; not_match ] );
  i = k([ 1; find(not_match)+1 ]);

end


% Output.
if( nargout>2 )

  % Compute and resort output index vector J so that A = B(J)
  % and also so that the duplicate entries in B points to the
  % correct entries in B.
  [i,jj] = sort( i );
  jjinv(jj) = 1:numel( jj );
  j = jjinv(j)';

elseif( nargout==1 )

  i = sort( i );
end

b = a( i, : );

if( dim==2 )
  b = b';
end
