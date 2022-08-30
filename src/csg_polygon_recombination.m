function [ polygons ] = csg_polygon_recombination( polygons, is_warn, tol )
%CSG_POLYGON_RECOMBINATION Recombine and tesselate polygons.
%
%   [ POLYGONS ] = CSG_POLYGON_RECOMBINATION( POLYGONS, IS_WARN, TOL )
%   Recombine co-linear/planar line segments and polygons. The IS_WARN
%   flag (default false) indicates is polygon recombination is skipped
%   due to differing identities.

% Initial version 180207.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help csg_polygon_recombination, return, end

if( nargin<3 )
  tol = 1e-5;
end
if( nargin<2 )
  is_warn = false;
end

ids = vertcat( polygons{1,:} );   % Polygon identities.
if( ~all(ids(1)==ids) )    % Only process if all identities are identical.
  if( is_warn )
    warning( 'polygon_recombination: no recombination attempted since identities differ.' )
  end
  return
end


dim = size(polygons{3,1},2);
n0  = size(polygons,2);
n   = n0;
is_changed = true;
n_cnt = 0;
while( n>1 && is_changed )
  is_changed = false;

  for i=1:n
    for j=1:n
      if( i~=j )

        if( dim==3 )
          joined_poly = l_join_polygons( polygons(:,i), polygons(:,j), tol );
        else
          joined_poly = l_join_line_segments( polygons(:,i), polygons(:,j) );
        end

        if( ~isempty(joined_poly) )
          polygons(:,[i,j]) = [];
          polygons = [ polygons, joined_poly ];
          n = size(polygons,2);
          is_changed = true;
          break;
        end

      end
    end
    if( is_changed ), break; end
  end

  n_cnt = n_cnt + 1;
  if( n_cnt>n0 )
    warning( 'recombine_polygons: max while loop counter reached.' )
    break;
  end
end


%------------------------------------------------------------------------------%
function [ joined_polygon ] = l_join_polygons( polygon_i, polygon_j, tol )

joined_polygon = [];


vi = polygon_i{4};
vj = polygon_j{4};
d  = polygon_i{6} + polygon_j{6};   % Mean diameter.
vi = vi/d;   % Scaling ~1.
vj = vj/d;
nr = polygon_i{3}';   % Reference normal.
% clf, plot(vi(:,1),vi(:,2),'rx'), hold on
% plot(vj(:,1),vj(:,2),'bo')

% Loop over vertices.
nvi = size(vi,1);
nvj = size(vj,1);
for iv1=1:nvi
  for jv1=1:nvj

    if( abs(vi(iv1,1)-vj(jv1,1))<tol && ...
        abs(vi(iv1,2)-vj(jv1,2))<tol && ...
        abs(vi(iv1,3)-vj(jv1,3))<tol )   % Matching points iv1==jv1.
      iv2 = mod(iv1,  nvi) + 1;
      jv0 = mod(jv1-2,nvj) + 1;
      jv2 = mod(jv1,  nvj) + 1;

      if( abs(vi(iv2,1)-vj(jv0,1))<tol && ...
          abs(vi(iv2,2)-vj(jv0,2))<tol && ...
          abs(vi(iv2,3)-vj(jv0,3))<tol )   % Matching points iv2==jv0.

        % Normal 1: iv0-iv1/jv1-jv2
        iv0 = mod(iv1-2,nvi) + 1;
        % p1 = [ vi([iv0,iv1],:); vj(jv2,:) ];
        n1 = cross( vi(iv1,:)-vi(iv0,:), vj(jv2,:)-vi(iv0,:) );

        % Normal 2: jvn-jv0/iv2-iv3
        jvn = mod(jv0-2,nvj) + 1;
        iv3 = mod(iv2,  nvi) + 1;
        % p2 = [ vj([jvn,jv0],:); vi(iv3,:) ];
        n2 = cross( vj(jv0,:)-vj(jvn,:), vi(iv3,:)-vj(jvn,:) );

        % Test for convexity.
        dnnr = [ n1*nr, n2*nr ];

        is_convex   = dnnr>tol;
        is_parallel = abs(dnnr)<=tol;

        if( all(is_convex | is_parallel) )

          % Points iv2 -> iv0 (skip iv1).
          ix = mod( iv2+[1:(nvi-1)]-1 -1, nvi ) + 1;
          if( is_parallel(2) )
            ix(1) = [];
          end

          % Points jv1 -> jvn (skip jv0).
          jx = mod( jv1+[1:(nvj-1)]-1 -1, nvj ) + 1;
          if( is_parallel(1) )
            jx(1) = [];
          end

          vij = [ polygon_i{4}(ix,:); polygon_j{4}(jx,:) ];
          if( size(vij,1)<3 )
            error( 'Could not form joined polygon.' )
          end
          joined_polygon = csg_op( vij, polygon_i{1}, 'b' );

          % plot_polygons( polygon_i{4}, polygon_j{4}, vij )
          return
        end

      elseif( abs(vi(iv2,1)-vj(jv2,1))<tol && ...
              abs(vi(iv2,2)-vj(jv2,2))<tol && ...
              abs(vi(iv2,3)-vj(jv2,3))<tol )   % Matching points iv2==jv2.

        is_identical = norm(vi(:)-vj(:)) < tol;
        if( is_identical )
          joined_polygon = polygon_i;
        else
          warning( 'Polygon orientation seems incorrect.' )
        end

      end
    end

  end
end


%------------------------------------------------------------------------------%
function plot_polygons( vi, vj, vij )

clf
hold on

plot3( vi( [1:end,1],1), vi( [1:end,1],2), vi( [1:end,1],3), 'g-',  'linewidth', 4 )
n1 = cross( vi(2,:)-vi(1,:), vi(3,:)-vi(1,:) );
c1 = mean(vi);
quiver3( c1(1), c1(2), c1(3), n1(1), n1(2), n1(3), 'color', 'g' )

plot3( vj( [1:end,1],1), vj( [1:end,1],2), vj( [1:end,1],3), 'r-',  'linewidth', 2 )
n2 = cross( vj(2,:)-vj(1,:), vj(3,:)-vj(1,:) );
c2 = mean(vj);
quiver3( c2(1), c2(2), c2(3), n2(1), n2(2), n2(3), 'color', 'r' )

if( nargin>=3 )
  plot3( vij([1:end,1],1), vij([1:end,1],2), vij([1:end,1],3), 'b--', 'linewidth', 1 )
  n3 = cross( vij(2,:)-vij(1,:), vij(3,:)-vij(1,:) );
  c3 = mean(vij);
  quiver3( c3(1), c3(2), c3(3), n3(1), n3(2), n3(3), 'color', 'b' )

  for k=1:size(vij,1)
    text( vij(k,1), vij(k,2), vij(k,3), num2str(k) )
  end
  vij
end

axis equal
view(3)
rotate3d('on')

pause

%------------------------------------------------------------------------------%
function [ joined_line ] = l_join_line_segments( line1, line2, TOL )

if( nargin<3 )
  TOL = sqrt(eps);
end

v1 = line1{4};
v2 = line2{4};
if( norm(v1(2,:)-v2(1,:))<TOL )

  v1(2,:) = v2(2,:);
  joined_line = csg_op( v1, line1{1}, 'b' );

elseif( norm(v2(2,:)-v1(1,:))<TOL )

  v2(2,:) = v1(2,:);
  joined_line = csg_op( v2, line1{1}, 'b' );

else
  joined_line = [];
end
