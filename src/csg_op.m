function [ c, d, stat ] = csg_op( a, b, op, do_merge, tol )
%CSG_OP Perform CSG operations.
%
%   [ C, (D), STAT ] = CSG_OP( A, B, OP, DO_MERGE, TOL ) Applies CSG
%   operation OP to trees A and B, where OP is a string character
%   indicating operation to perform ('+' union/join, '-' subtract, '&'
%   intersect). If DO_MERGE is true (default) the resulting CSG trees
%   are merged, pruned, and returned in C, or if two outputs are
%   requested the merged CSG trees from A and B are returned in C and
%   D, respectively. The status return code STAT is zero if the
%   operation has been applied sucessfully, and a positive integer to
%   indicate error in which case the original inputs are returned
%   unmodified.
%
%   With OP equal to 'b' then the CSG build operation will be called
%   on the polygons or vertices in A (with optional polygon identity
%   argument in B). OP 'e' extracts and returns all polygons from the
%   CSG tree A (and corresponding CSG node normals in D). OP 'u'
%   updates the polygon vertices with those in B, assuming same order
%   as the extract operation. If only two input arguments are given,
%   then B will be used as OP. OP 'i' will add the contents of B to
%   the polygon identities. OP 'p' will prune the CSG tree A, that is
%   remove empty CSG nodes. OP 'v' and 'w' plots and visualizes the
%   polygons and CSG tree nodes in A, respectively.
%
%   TOL (default 1e-7) specifies the tolerance for determining point
%   relative to splitting plane positions:
%
%       t =  dot(plane_normal,point) - dot(plane_normal,point_in_plane)
%
%   so that t>tol is in front of the plane, and t<-tol is behind.
%
%   The CSG tree data structure are nested 1 x 5 cell arrays with the
%   following contents for each node:
%
%   { single point coordinates in node plane, plane normal,
%     nested front nodes, nested back nodes, polygons }
%
%   The polygon data structure is a 6 x n_polygons cell array where
%   the respective rows correspond to: polygon indices, number of
%   vertices per cell, normal, vertex coordinates, bounding
%   circle/sphere center, and radius.

% Initial version 171228.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help csg_op, return, end

if( nargin<5 )
  tol = 1e-7;
if( nargin<4 )
  do_merge = true;
if( nargin<3 )
  op = b;
end,end,end
stat = 0;


NRMAX = 500;
if( exist( 'OCTAVE_VERSION', 'builtin' ) )
  max_recursion_depth( NRMAX );
else
  set( 0, 'RecursionLimit', NRMAX );
end


% Apply CSG operation.
switch( op )

  case {'+'}          % Join, A + B
    [c,d,stat] = union( a, b, do_merge, tol );

  case {'-'}          % Subtract, A - B
    [c,d,stat] = subtract( a, b, do_merge, tol );

  case {'&'}          % Intersect, A & B
    [c,d,stat] = intersect( a, b, do_merge, tol );

  case {'b','B'}      % Build

    % Build CSG tree from polygons A
    if( iscell(a) && size(a,1)==6 )
      c = node( a, tol );

    % Build polygon from vertices A and identity argument B
    elseif( isnumeric(a) && size(a,1)>=2 && any(size(a,2)==[2,3]) )
      c = build_polygon( a, b );
    else
      error( 'A input not supported for build operation.' )
    end

  case {'p','P'}      % Prune CSG tree
    c = prune( a, tol );

  case {'i','I','id','ID'}    % Add B to polygon identities
    c = add_to_polygon_id( a, b );

  case {'e','E'}      % Extract polygons (and normals)
    [c,d] = extract_polygons( a );

  case {'u','U'}      % Update polygon vertices
    c = update_polygons( a, b );

  case {'v','V'}      % Plot and visualize polygons
    if( nargin<3 )
      b = 1;
    end
    plot_polygons( a, b );

  case {'w','W'}      % Plot and visualize CSG nodes
    plot_csg( a );

  case {'wp','WP'}    % Plot and visualize polygons of CSG nodes
    plot_csgp( a );

  otherwise
    error( [ char(op), ' is not a valid operation.'] )
end


if( stat~=0 )
  c = a;
  d = b;
end


%------------------------------------------------------------------------------%
% CSG Tree operations.
%------------------------------------------------------------------------------%
function [ c, d, stat ] = union( a, b, do_merge, tol )
% Union of A and B (A + B).

METHOD = 1;

a = clip( a, b, false, tol );

if( METHOD==1 )
  b = clip( b, a, true, tol );
else
  b = clip( b, a, false, tol );
  b = invert( b );   % Remove coplanar polygons/colinear lines from b.
  b = clip( b, a, false, tol );
  b = invert( b );
end


% Merge CSG trees.
if( do_merge )
  c = merge( a, b );
  c = prune( c, tol );
  d = [];
else
  c = a;
  d = b;
end

stat = 0;

%------------------------------------------------------------------------------%
function [ c, d, stat ] = subtract( a, b, do_merge, tol )
% Subtraction of B from A (A - B).

METHOD = 2;

[tmp,n_nodes_a0] = count_polygons( a );

a = invert( a );
a = clip( a, b, false, tol );

if( METHOD==1 )
  b = clip( b, a, true, tol );
else
  b = clip( b, a, false, tol );
  b = invert( b );
  b = clip( b, a, false, tol );
  b = invert( b );
end


[tmp,n_nodes_a] = count_polygons( a );
n_polygons_b    = count_polygons( b );
if( n_polygons_b==0 && n_nodes_a==n_nodes_a0 )
  stat = 1;
  c = [];
  d = [];
  return;
end


% Merge CSG trees.
if( do_merge )
  c = build( a, extract_polygons(b) );
  c = invert( c );
  c = prune( c, tol );
  d = [];
  if( isempty(c) )
    stat = 1;
    return;
  end
else
  c = invert( a );
  d = invert( b );
end
stat = 0;

%------------------------------------------------------------------------------%
function [ c, d, stat ] = intersect( a, b, do_merge, tol )
% Intersection of A and B (A & B).

a = invert( a );
b = clip( b, a, false, tol );
b = invert( b );
a = clip( a, b, false, tol );
b = clip( b, a, false, tol );  % FIXME: Unnecessary?


n_polygons_a = count_polygons( a );
n_polygons_b = count_polygons( b );
if( (n_polygons_a + n_polygons_b)==0 )
  stat = 1;
  c = [];
  d = [];
  return;
end


% Merge CSG trees.
if( do_merge )
  c = merge( a, b, n_polygons_a, n_polygons_b );
  c = invert( c );
  c = prune( c, tol );
  d = [];
else
  c = invert( a );
  d = invert( b );
end
stat = 0;


%------------------------------------------------------------------------------%
function [ a ] = clip( a, b, rem_cfront, tol )
% Remove all polygons in A that are inside B (clip B from A, A [-] B).

if( nargin<4 )
  tol = 1e-5;
if( nargin<3 )
  rem_cfront = false;
end,end

polygons = a{5};
if( ~isempty(polygons) )

  polygons = clip_polygons( polygons, b, rem_cfront, tol );

  % Re-combine/tesellate split polygons.
  n_polygons = size(polygons,2);
  if( n_polygons>1 )
    polygons = csg_polygon_recombination( polygons );
  end

  a{5} = polygons;
end

% Process node to the front.
if( ~isempty(a{3}) )
  a{3} = clip( a{3}, b, rem_cfront, tol );
end

% Process node to the back.
if( ~isempty(a{4}) )
  a{4} = clip( a{4}, b, rem_cfront, tol );
end

%------------------------------------------------------------------------------%
function [ a ] = invert( a )
% Invert A (converts solid space to empty space and vice versa).

% Reorient plane normal.
if( ~isempty(a{2}) )
  a{2} = -a{2};
end

% Exchange front and back nodes.
front = a{3};
a{3}  = a{4};
a{4}  = front;

% Reorient polygon vertices and normal.
polygons = a{5};
if( ~isempty(polygons) )
  polygons = invert_polygons( polygons );

  a{5} = polygons;
end

% Invert front and back nodes (recursion).
if( ~isempty(a{3}) )
  a{3} = invert( a{3} );
end
if( ~isempty(a{4}) )
  a{4} = invert( a{4} );
end

%------------------------------------------------------------------------------%
function [ c ] = merge( a, b, n_polygons_a, n_polygons_b )
% Merge CSG trees A and B.

if( nargin<4 )
  n_polygons_b = count_polygons( b );
if( nargin<3 || isempty(n_polygons_a) )
  n_polygons_a = count_polygons( a );
end,end

if( n_polygons_a==0 )
  c = b;
elseif( n_polygons_b==0 )
  c = a;
elseif( n_polygons_a < n_polygons_b )
  c = build( a, extract_polygons(b) );
else
  c = build( b, extract_polygons(a) );
end

%------------------------------------------------------------------------------%
function [ a ] = prune( a, tol )
% Prune and remove empty nodes from CSG tree A.

if( ~isempty(a{3}) )
  a{3} = prune( a{3}, tol );
end
if( ~isempty(a{4}) )
  a{4} = prune( a{4}, tol );
end

if( isempty(a{5}) )

  if( isempty(a{3}) && isempty(a{4}) )
    a = [];
  elseif( ~isempty(a{3}) &&  isempty(a{4}) )
    a = a{3};
  elseif(  isempty(a{3}) && ~isempty(a{4}) )
    a = a{4};
  else
    a = build( a{3}, extract_polygons(a{4}), tol );
  end
end

%------------------------------------------------------------------------------%
function [ a ] = node( polygons, tol )
% Create an empty CSG node (and populate with polygons).

if( nargin<2 )
  tol = 1e-5;
end

a = cell(1,5);   % a = { plane point, plane normal, front , back, polygons };
a{5} = {};

if( nargin>=1 )
  a = build( a, polygons, tol );
end

%------------------------------------------------------------------------------%
function [ a ] = build( a, polygons, tol )
% Build a BSP CSG tree out of 3D polygons or 2D line segments. When
% called on an existing tree A, the new polygons are filtered down to
% the bottom of the tree and become new nodes there. Each set of
% polygons is partitioned using the first polygon (no heuristic is
% used to pick a good split).

if( nargin<3 )
  tol = 1e-5;
if( nargin<2 || isempty(polygons) )
  return
end,end

if( isempty(a{1}) )   % Default splitting plane.
  [a{1},a{2}] = get_plane( polygons );
end
point  = a{1};
normal = a{2};


n_polygons = size(polygons,2);
if( n_polygons>1 )

  [front,back,ix_split] = ...
      sort_polygons( polygons, point, normal, true, tol );

else

  front = {};
  back  = {};
  ix_split = 1;

end

if( ~isempty(ix_split) )
  for i=ix_split

    [front,back,cfront,cback] = ...
        split_polygon( polygons(:,i), point, normal, front, back, tol );

    if( ~isempty(cfront) )
      a{5} = [ a{5}, cfront ];
    end
    if( ~isempty(cback) )
      a{5} = [ a{5}, cback ];
    end

  end
end


if( ~isempty(front) )
  if( isempty(a{3}) )
    a{3} = node();
  end
  a{3} = build( a{3}, front, tol );
end

if( ~isempty(back) )
  if( isempty(a{4}) )
    a{4} = node();
  end
  a{4} = build( a{4}, back, tol );
end



%------------------------------------------------------------------------------%
% CSG Tree and Polygon operations.
%------------------------------------------------------------------------------%
function [ polygons ] = clip_polygons( polygons, b, rem_cfront, tol )
% Recursively remove all polygons/segments inside B.
%
% Loop over all polygons, and split into front and back parts
% according to all polygons in B.

if( isempty(b{1}) || isempty(polygons) )
  return
end
if( nargin<4 )
  tol = 1e-5;
if( nargin<3 )
  rem_cfront = false;
end,end
point  = b{1};
normal = b{2};


% Inlined and unrolled sort_polygons_by_bound call.
n_polygons = size(polygons,2);
front = {};
back  = {};
if( n_polygons==2 )

  normal = normal';
  pdn    = point * normal;
  dists  = [ polygons{5,1} * normal - pdn , ...
             polygons{5,2} * normal - pdn ];
  normal = normal';

  radii = [ polygons{6,1}, polygons{6,2} ];
  ix_front = dists >  radii;
  ix_back  = dists < -radii;

  ix_split = find( ~(ix_front | ix_back) );
  ix_front = find( ix_front );
  if( ~isempty(ix_front) )
    front = polygons(:,ix_front);
  end
  ix_back = find( ix_back );
  if( ~isempty(ix_back) )
    back = polygons(:,ix_back);
  end

elseif( n_polygons==1 )

  center = polygons{5};
  radii  = polygons{6};
  dist   = normal * (center - point)';
  if( dist > radii )
    front = polygons;
    ix_split = [];
  elseif( dist < -radii )
    back = polygons;
    ix_split = [];
  else
    ix_split = 1;
  end

elseif( n_polygons==3 )

  normal = normal';
  pdn    = point * normal;
  dists  = [ polygons{5,1} * normal - pdn , ...
             polygons{5,2} * normal - pdn , ...
             polygons{5,3} * normal - pdn ];
  normal = normal';

  radii = [ polygons{6,1}, polygons{6,2}, polygons{6,3} ];
  ix_front = dists >  radii;
  ix_back  = dists < -radii;

  ix_split = find( ~(ix_front | ix_back) );
  ix_front = find( ix_front );
  if( ~isempty(ix_front) )
    front = polygons(:,ix_front);
  end
  ix_back = find( ix_back );
  if( ~isempty(ix_back) )
    back = polygons(:,ix_back);
  end

else   % n_polygons>3

  [ix_front,ix_back,ix_split] = ...
      sort_polygons_by_bound( polygons, point, normal, tol );

  front = polygons(:,ix_front);
  back  = polygons(:,ix_back);
  ix_split = find( ix_split );

end


if( ~isempty(ix_split) )
  for i=ix_split

    [front,back,cfront,cback] = ...
        split_polygon( polygons(:,i), point, normal, front, back, tol );

    if( ~isempty(cfront) )
      if( rem_cfront )   % Remove coplanar front nodes.
        back = [ back, cfront ];
      else
        front = [ front, cfront ];
      end
    end
    if( ~isempty(cback) )
      back = [ back, cback ];
    end
  end
end


if( ~isempty(b{3}) )
  front = clip_polygons( front, b{3}, rem_cfront, tol );
end
if( ~isempty(b{4}) )
  back = clip_polygons( back, b{4}, rem_cfront, tol );
else
  back = {};
end


if( ~isempty(front) )
  polygons = front;
  if( ~isempty(back) )
    polygons = [ polygons, back ];
  end
else
  if( ~isempty(back) )
    polygons = back;
  else   % Delete nodes behind plane of b.
    polygons = cell(6,0);
  end
end

%------------------------------------------------------------------------------%
function [ polygons, normals ] = extract_polygons( a )
% Extracts all polygons/segments and corresponding normals from CSG tree A.

polygons = a{5};
normals  = repmat(a{2},size(polygons,2),1);

if( ~isempty(a{3}) )   % Front.
  [front_polygons,front_normals] = extract_polygons( a{3} );
  if( isempty(polygons) )
    polygons = front_polygons;
  else
    polygons = [ polygons, front_polygons ];
  end
  normals = [ normals; front_normals ];
end

if( ~isempty(a{4}) )   % Back.
  [back_polygons,back_normals] = extract_polygons( a{4} );
  if( isempty(polygons) )
    polygons = back_polygons;
  else
    polygons = [ polygons, back_polygons ];
  end
  normals = [ normals; back_normals ];
end

%------------------------------------------------------------------------------%
function [ a, v ] = update_polygons( a, v )
% Updates polygons vertices in CSG tree A. Inserts vertices
% from V in the same order as for polygon extraction.

if( ~isempty(a{5}) )
  for i=1:size(a{5},2)
    n = size(a{5}{4,i},1);
    a{5}{4,i}  = v(1:n,:);
    v(1:n,:) = [];
  end
end

if( ~isempty(a{3}) )   % Front.
  [a{3},v] = update_polygons( a{3}, v );
end

if( ~isempty(a{4}) )   % Back.
  [a{4},v] = update_polygons( a{4}, v );
end

%------------------------------------------------------------------------------%
function [ n_polygons, n_nodes ] = count_polygons( a )
% Count the total number of polygons in CSG tree A.

n_polygons = size(a{5},2);
n_nodes    = n_polygons > 0;

if( ~isempty(a{3}) )   % Front.
  [n_front_polygons,n_front_nodes] = count_polygons( a{3} );

  n_polygons = n_polygons + n_front_polygons;
  n_nodes    = n_nodes + n_front_nodes;
end

if( ~isempty(a{4}) )   % Back.
  [n_back_polygons,n_back_nodes] = count_polygons( a{4} );

  n_polygons = n_polygons + n_back_nodes;
  n_nodes    = n_nodes + n_back_nodes;
end

%------------------------------------------------------------------------------%
function [ a ] = add_to_polygon_id( a, id )
% Add ID to the polygons identity property.

polygons = a{5};
if( ~isempty(polygons) )
  n_polygons = size(polygons,2);
  for i=1:n_polygons
    id_i = polygons{1,i}(:,end);
    polygons{1,i}(:,end) = id_i + id;
  end
  a{5} = polygons;
end

% Process node to the front.
if( ~isempty(a{3}) )
  a{3} = add_to_polygon_id( a{3}, id );
end

% Process node to the back.
if( ~isempty(a{4}) )
  a{4} = add_to_polygon_id( a{4}, id );
end


%------------------------------------------------------------------------------%
% Polygon operations.
%------------------------------------------------------------------------------%
function [ front, back, cfront, cback ] = split_polygon( polygon, a_point, a_normal, front, back, tol )
% Split POLYGON by the plane of A (if needed), then put the entire polygon
% or polygon fragments in the appropriate lists. If the coplanar output
% argument not requested, then coplanar polygons go into either front or
% back depending on their orientation with respect to the plane of A.

if( nargin<6 )
  tol = 1e-5;
end

% Classification.
COPLANAR = 0;
FRONT    = 1;
BACK     = 2;
SPANNING = 3;

p_normal   = polygon{3};
p_vertices = polygon{4};
n_v = size(p_vertices,1);
dim = size(p_vertices,2);

ndp = a_normal * a_point';
t = sum( repmat(a_normal,n_v,1) .* p_vertices, 2 ) - ndp;
types = double(t>tol) + BACK*double(t<-tol);   % FRONT & BACK
polygon_type = max( types );
if( any(types==FRONT) && any(types==BACK) )
  polygon_type = SPANNING;
end

% Put the polygon in the correct list, splitting it when necessary.
cfront = {};
cback  = {};
switch( polygon_type )

  case COPLANAR
    if( nargout>2 )
      if( (a_normal*p_normal')>0 )
        cfront = [ cfront, polygon ];
      else
        cback = [ cback, polygon ];
      end
    else
      if( (a_normal*p_normal')>0 )
        front = [ front, polygon ];
      else
        back = [ back, polygon ];
      end
    end

  case FRONT
    front = [ front, polygon ];

  case BACK
    back = [ back, polygon ];

  case SPANNING
    f = [];
    b = [];
    for i=1:n_v
      j  = mod(i,n_v) + 1;
      ti = types(i);
      tj = types(j);
      vi = p_vertices(i,:);
      vj = p_vertices(j,:);

      if( ti ~= BACK )
        f = [ f; vi ];
      end
      if( ti ~= FRONT )
        b = [ b; vi ];
      end
      if( (ti+tj)==SPANNING )
        vji = ( vj - vi );
        t = ( ndp - a_normal*vi' ) / ( a_normal*vji' );
        v = vi + t*vji;
        f = [ f; v ];
        b = [ b; v ];
      end
    end

    if( size(f,1)>=dim )
      front = [ front, build_polygon( f, polygon{1} ) ];
    end
    if( size(b,1)>=dim )
      back = [ back, build_polygon( b, polygon{1} ) ];
    end

end

%------------------------------------------------------------------------------%
function [ polygons ] = invert_polygons( polygons )
% Invert POLYGONS by reversing normal and vertex order.

n_polygons = size(polygons,2);
for i=1:n_polygons

  % Reverse normal.
  polygons{3,i} = -polygons{3,i};

  % Reverse vertices.
  polygons{4,i} = polygons{4,i}(end:-1:1,:);
end


%------------------------------------------------------------------------------%
function [ polygon ] = build_polygon( vertices, identity )
% Construct polygon/line segment from VERTICES.

if( nargin<2 )
  identity = [];
end

% Polygon normal vector.
p1 = vertices(1,:);
p2 = vertices(2,:);
if( size(vertices,2)==2 )   % 2D - Line segment.

  vertices = vertices(1:2,:);
  d  = vertices(2,:) - vertices(1,:);   % Line direction/tangent vector.
  n  = [ d(2), -d(1) ];   % Normal vector pointing to the right of the line direction.

else   % 3D - Polygon

  p3 = vertices(3,:);
  n  = cross( p2-p1, p3-p1 );

end
normal = n./sqrt(sum(n.^2));

% Center and radius of bounding circle/sphere.
center = ( max(vertices) + min(vertices) )/2;   % Polygon (bounding box) center.
TOL    = eps*1e-3;
radius = sqrt(max( sum( (vertices - ...
              repmat(center,size(vertices,1),1)).^2, 2 ) )) + TOL;


% Define polygon
n_v_p = size(vertices,1);   % Number of vertices per polygon.
polygon = { identity ;
            n_v_p    ;
            normal   ;
            vertices ;
            center   ;
            radius   };

%------------------------------------------------------------------------------%
function [front,back,ix_split] = sort_polygons( polygons, point, normal, use_bound, tol )
% Sort POLYGONS to FRONT and BACK of the plane with POINT and NORMAL.
% The output argument IX_SPLIT is an index vector indicating polygons
% that are split by the plane (must be done in a separate routine).
% The USE_BOUND flag specified to also use bounding sphere (default)
% to speed up front/back testing.

if( nargin<5 )
  tol = 1e-5;
end
if( nargin<4 || use_bound )
  [ix_front,ix_back,ix_split] = sort_polygons_by_bound( polygons, point, normal, tol );
else
  n = size(polygons,2);
  ix_front = zeros(1,n);
  ix_back  = ix_front;
  ix_split = ones(1,n);
end

% Classification.
COPLANAR = 0;
FRONT    = 1;
BACK     = 2;
SPANNING = 3;

ix_split = find(ix_split);
n_v_p = [ polygons{2,ix_split} ];   % Number of vertices for each polygon.
v_all = vertcat( polygons{4,ix_split} );
is = cumsum( [0,n_v_p(1:end-1)] );

polygon_type = ix_front + BACK*ix_back;
while( any(n_v_p) )
  n_v_i = n_v_p(find(n_v_p,1,'first'));
  ix = find(n_v_p==n_v_i);

  ixv = repmat(is(ix),n_v_i,1) + repmat([1:n_v_i]',1,length(ix));
  v = v_all(ixv,:);
  n_v_tot = size(v,1);

  t = reshape( sum( repmat(normal,n_v_tot,1) .* ...
                    ( v - repmat(point,n_v_tot,1) ), 2 ), n_v_i, n_v_tot/n_v_i );

  types = double( t > tol ) + BACK*double( t < -tol );   % FRONT & BACK
  polygon_type(ix_split(ix)) = max( types );
  msk = any(types==FRONT) & any(types==BACK);
  polygon_type(ix_split(ix(msk))) = SPANNING;

  n_v_p(ix) = 0;
end

front = polygons( :, polygon_type==FRONT );
back  = polygons( :, polygon_type==BACK  );
ix_split = find(polygon_type==SPANNING | polygon_type==COPLANAR);

%------------------------------------------------------------------------------%
function [ix_front,ix_back,ix_split] = sort_polygons_by_bound( polygons, point, normal, tol )
% Sort POLYGONS to front and back of the plane with POINT and NORMAL
% by using the polygons bounding sphere/circle. Returns the index
% logical masks IX_FRONT, IX_BACK, and IX_SPLIT.

if( nargin<4 )
  tol = 1e-5;
end

n = size(polygons,2);
centers = vertcat( polygons{5,:} );
radii   = [ polygons{6,:} ];

dists    = sum( repmat(normal,n,1) .* ...
                (centers - repmat(point,n,1)), 2 )';

ix_front = dists >  radii + tol;
ix_back  = dists < -radii - tol;
ix_split = ~(ix_front | ix_back);

%------------------------------------------------------------------------------%
function [ point, normal ] = get_plane( polygons, idx )
% Get point and normal from polygon IDX (or last polygon).

if( nargin<2 )
  idx = 1;
end

n_polygons = size(polygons,2);
i = min(idx,n_polygons);

vertices = polygons{4,i};
point    = vertices(1,:);
normal   = polygons{3,i};



%------------------------------------------------------------------------------%
% Utility and help functions
%------------------------------------------------------------------------------%
function plot_csg( a, color )
% Plot CSG tree A.

if( nargin<2 )
  color = 0.9*[1,1,1];
end
LW = 2;

point  = a{1};
normal = a{2};

hold on
s = 1/4;
quiver( point(1), point(2), normal(1), normal(2), s, 'linewidth', LW )
plot( point(1), point(2), 'ko' )

if( ~isempty(a{3}) )   % Front.
  plot_csg( a{3}, color )
end

if( ~isempty(a{4}) )   % Back.
  plot_csg( a{4}, color )
end

%------------------------------------------------------------------------------%
function plot_csgp( a, color )
% Plot polygons of CSG nodes in A (pause after each node).

if( nargin<2 )
  color = 0.9*[1,1,1];
end
LW = 2;

if( ~isempty(a{5}) )
  b    = a;
  b{3} = [];
  b{4} = [];
  hp = plot_polygons( b, 1, 'r' );
  pause
  set( hp, 'facecolor', 0.9*[1,1,1] )
end

if( ~isempty(a{3}) )   % Front.
  plot_csgp( a{3}, color )
end

if( ~isempty(a{4}) )   % Back.
  plot_csgp( a{4}, color )
end

%------------------------------------------------------------------------------%
function [ hp ] = plot_polygons( a, style, color )
% Plot polygons in CSG tree A.

if( nargin<3 )
  color = 0.9*[1,1,1];
end
LW = 2;

polygons = extract_polygons( a );
n_polygons = size(polygons,2);
n_sdim = size(polygons{3,1},2);

hold on
hp = [];
for i=1:n_polygons
  v = polygons{4,i};
  pc = mean(v);
  if( n_sdim==2 )
    x = v(1,1);
    y = v(1,2);
    u = v(2,1)-x;
    v = v(2,2)-y;
    quiver( x, y, u, v, 1, 'linewidth', LW, 'maxheadsize', 10 )
    text( pc(1), pc(2), num2str(i) )
  else
    f = 1:size(v,1);
    hp = [ hp, patch( 'faces', f, 'vertices', v, 'facecolor', color, 'linewidth', LW ) ];
    text( pc(1), pc(2), pc(3), num2str(i) )
    if( style==2 )
      set( hp, 'facecolor', 'g' )
      h = plot3( v(:,1), v(:,2), v(:,3), 'r.', 'markersize', 20 );
    end
  end

  if( style==2 )
    axis equal
    xlabel( 'x ')
    ylabel( 'y ')
    if( n_sdim==3 )
      zlabel( 'z ')
      alpha( 0.5 )
      view(3)
      rotate3d('on')
    else
      grid on
    end
    pause
    if( n_sdim==3 )
      set( hp, 'facecolor', color )
    end
    delete( h )
  end
end
axis equal
xlabel( 'x ')
ylabel( 'y ')
if( n_sdim==3 )
  zlabel( 'z ')
  alpha( 0.5 )
  view(3)
  rotate3d('on')
else
  grid on
end

if( ~nargout )
  clear hp
end
