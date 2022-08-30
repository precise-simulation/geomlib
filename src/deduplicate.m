function [ b, i, j ] = deduplicate( a, dim, tol, mode )
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


%% UNIT TESTS

%!assert( isequal( deduplicate([0 0;0 1;0 0;0 1]), [0 0;0 1] ) );
%!assert( isequal( deduplicate([0 0;0 1;0 0;0 1;1 0],1), [0 0;0 1;1 0] ) );
%!assert( isequal( deduplicate([0 0;0 1;0 0;0 1]',2), [0 0;0 1]' ) );
%!assert( isequal( deduplicate([0 0;0 1+1e4*eps;0 0;0 1]), [0 0;0 1+1e4*eps;0 1] ) );
%!assert( isequal( deduplicate([0 0;0 1+1e4*eps;0 0;0 1],1,3e4*eps), [0 0;0 1+1e4*eps] ) );

%!assert( isequal( deduplicate([0 0;0 1;0 0;0 1],1,eps*1e3,0), [0 0;0 1] ) );
%!assert( isequal( deduplicate([0 0;0 1;0 0;0 1;1 0],1,eps*1e3,0), [0 0;0 1;1 0] ) );
%!assert( isequal( deduplicate([0 0;0 1;0 0;0 1]',2,eps*1e3,0), [0 0;0 1]' ) );
%!assert( isequal( deduplicate([0 0;0 1+1e4*eps;0 0;0 1],1,eps*1e3,0), [0 0;0 1+1e4*eps;0 1] ) );
%!assert( isequal( deduplicate([0 0;0 1+1e4*eps;0 0;0 1],1,3e4*eps,0), [0 0;0 1+1e4*eps] ) );

%!test % Test 1
%! a = [0 0;0.5 0.5;1 1;0 0;0 -1;-1 -1;1 1;0 -1;-1 -1;1 1;-1 -1;0 -1;-2 -2];
%! b_ref = [0 0;0.5 0.5;1 1;0 -1;-1 -1;-2 -2];
%!
%! [b,i,j] = deduplicate( a );
%! assert( isequal( b, b_ref ) )
%! assert( isequal( b, a(i,:) ) )
%! assert( isequal( a, b(j,:) ) )
%!
%! [b,i,j] = deduplicate( a, 1, eps, 0 );
%! assert( isequal( b, b_ref ) )
%! assert( isequal( b, a(i,:) ) )
%! assert( isequal( a, b(j,:) ) )

%!test % Test 2
%! a = [0 0;0.5 0.5;1 1;0 eps*1e2;0 -1;-1 -1;1 1;0 -1-eps*2e2;-1 -1;1 1;1.1 1.2;-1.3 -1.4];
%! b_ref = [0 0;0.5 0.5;1 1;0 -1;-1 -1;1.1 1.2;-1.3 -1.4];
%!
%! [b,i,j] = deduplicate( a );
%! assert( isequal( b, b_ref ) )
%! assert( isequal( b, a(i,:) ) )
%! assert( all(all(abs( a-b(j,:) )<eps*1e3)) )
%!
%! [b,i,j] = deduplicate( a, 1, eps*1e3, 0 );
%! assert( isequal( b, b_ref ) )
%! assert( isequal( b, a(i,:) ) )
%! assert( all(all(abs( a-b(j,:) )<eps*1e3)) )
