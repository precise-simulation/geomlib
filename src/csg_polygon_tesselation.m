function [ vertices, faces, ids ] = csg_polygon_tesselation( polygons, tol )
%CSG_POLYGON_TESSELATION Recombine and tesselate polygons.
%
%   [ VERTICES, FACES, IDS ] = CSG_POLYGON_TESSELATION( POLYGONS, TOL )
%   Recombine and tesselate POLYGONS.

% Initial version 180319.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help csg_polygon_tesselation, return, end

if( nargin<2 )
  tol = 1e-5;
end

IS_DBG = false;

vertices = [];
faces    = [];
ids      = [];
if( isempty(polygons) || isempty(polygons{4}) )
  return;
end


% Group polygons.
ids = vertcat( polygons{1,:} );   % Polygon identities.
if( ~isempty(ids) )
  if( size(ids,2)==1 )
    [tmp,ind] = uunique(ids);
  else
    [tmp,ind] = uunique(ids(:,end-1:end));
  end
  if( size(ids,1)==1 )
    ind = 1;
  end
else
  ids = ones(size(polygons,2),1);
  ind = ids;
end

if( IS_DBG )
  clf
  subplot(1,2,1)
  plot_polygons( polygons )
end

ids_match = [];
for i=1:max(ind)
  [vertices_i,faces_i] = l_recombine_polygons( polygons(:,ind==i), tol );

  offset    = size(vertices,1);
  faces     = [ faces; faces_i+offset ];
  vertices  = [ vertices; vertices_i  ];
  ids_match = [ ids_match; repmat(ids(find(ind==i,1),:),size(faces_i,1),1) ];
end
ids = ids_match;

if( IS_DBG )
  faces
  subplot(1,2,2)
  patch( 'vertices', vertices, 'faces', faces, 'facecolor', 'g' )
  hold on
  n  = l_polygon_normal_area( vertices, faces );
  pc = [ mean(reshape(vertices(faces,1),size(faces)),2), ...
         mean(reshape(vertices(faces,2),size(faces)),2), ...
         mean(reshape(vertices(faces,3),size(faces)),2) ];
  for i=1:size(n,1)
    quiver3(pc(i,1),pc(i,2),pc(i,3),n(i,1),n(i,2),n(i,3))
  end
  view(3), pause
end


%------------------------------------------------------------------------------%
function [ vertices, t ] = l_recombine_polygons( polygons, tol )

A_TOL = 1e-9;

n_polygons = size(polygons,2);
if( n_polygons==1 )
  [vertices,t] = l_recombine_polygon( polygons, tol );
  return;
end

% Collect vertices and faces from all polygons.
vertices = [];
c_faces  = cell(1,n_polygons);
for i=1:n_polygons
  v_i = polygons{4,i};
  offset     = size(vertices,1);
  vertices   = [ vertices; v_i ];
  c_faces{i} = [1:size(v_i,1)] + offset;
end

% Deduplicate vertices.
[vertices,tmp,ind_v_orig] = deduplicate( vertices, 1, tol );
for i=1:length(c_faces)
  c_faces{i} = reshape( ind_v_orig(c_faces{i}), size(c_faces{i}) );
end

if( n_polygons==2 && size(vertices,1)==4 && ...
    length(c_faces{1})==3 && length(c_faces{2})==3 )

  t = vertcat( c_faces{:} );

  e_ext = l_get_external_edges( t );
  iv_rem = l_collinear_edge_vertices( vertices, e_ext );

  if( ~isempty(iv_rem) )
    if( length(iv_rem)>1 )
      vertices = [];
      t = [];
    end

    n1 = l_polygon_normal_area( vertices, t );
    vertices(iv_rem,:) = [];
    t  = 1:3;
    n2 = l_polygon_normal_area( vertices, t );
    if( dot(mean(n1),n2)<0 )
      t = 3:-1:1;
    end

  end

  return;
end

% Project vertices to 2D plane.
[v2,n1] = l_project_boundary( vertices, c_faces );

% Compute triangulation.
t = l_delaunay( v2 );

% Remove triangles with center outside original faces.
pc = [ mean(reshape(v2(t,1),size(t)),2), ...
       mean(reshape(v2(t,2),size(t)),2) ];
ix_t_keep = l_in_polygon( pc, v2, c_faces );

% Remove triangles with minimal area (<A_TOL*max area).
d12 = v2(t(:,2),:) - v2(t(:,1),:);
d13 = v2(t(:,3),:) - v2(t(:,1),:);
a = ( d12(:,1).*d13(:,2) - d12(:,2).*d13(:,1) )/2;
t(a<0,:) = t(a<0,3:-1:1);
ix_t_keep = ix_t_keep & abs(a)>max(abs(a))*A_TOL;
t = t(ix_t_keep,:);
% clf, patch( 'vertices', v2, 'faces', t, 'facecolor', [.9 .9 .9] ), pause

% Remove internal and collinear edge vertices.
e_ext = l_get_external_edges( t );
ind_v_keep = setdiff( unique(e_ext(:).'), l_collinear_edge_vertices( vertices, e_ext ) );

% Retriangulate.
v22 = v2(ind_v_keep,:);
t = l_delaunay( v22 );
% clf, patch( 'vertices', v22, 'faces', t, 'facecolor', [.9 .9 .9] ), pause

% Remove triangles with center outside original faces.
if( ~all(ix_t_keep) )
  pc = [ mean(reshape(v22(t,1),size(t)),2), ...
         mean(reshape(v22(t,2),size(t)),2) ];
  ix_t_keep = l_in_polygon( pc, v2, c_faces );
  t = t(ix_t_keep,:);
end

% Remove triangles with minimal area (<A_TOL*max area).
d12 = v22(t(:,2),:) - v22(t(:,1),:);
d13 = v22(t(:,3),:) - v22(t(:,1),:);
a = ( d12(:,1).*d13(:,2) - d12(:,2).*d13(:,1) )/2;
t(a<0,:) = t(a<0,3:-1:1);
t = t(abs(a)>max(abs(a))*A_TOL,:);

vertices = vertices(ind_v_keep,:);

% Reorient
n2 = l_polygon_normal_area( vertices, t );
ix_reorient = dot(repmat(n1,size(n2,1),1),n2,2)<0;
t(ix_reorient,:) = t(ix_reorient,3:-1:1);

%------------------------------------------------------------------------------%
function [ vertices_2d, nm ] = l_project_boundary( vertices, faces )

% Projection plane normal and point.
[n,a] = l_polygon_normal_area( vertices, faces );
n0 = n;
if( size(n,1)>1 )
  n = mean( n.*repmat(a,1,3) );
end
if( all(abs(n)<eps*1e1) || any(isnan(n)) )
  [tmp,i] = max(a);
  n = n0(i,:);
end
n = n/norm(n);
nm = mean(n,1);
p_c = mean( vertices );

% Project points to plane.
% for i=1:size(vertices,1)
%   p_i = vertices(i,:);
%   p_i = p_i - dot(n,(p_i-p_c))*n;
%   vertices(i,:) = p_i;
% end

% Orthonormal projection to 2D plane.
v = zeros(1,3);
i = 1;
while( all(abs(v)<eps*1e1) )
  if( i>size(vertices,1) )
    error( 'Failed to determine projection plane.' )
  end

  p_t = vertices(i,:);
  p_t = p_t - dot(n,(p_t-p_c))*n;   % Point 1 in plane.
  t = p_t - p_c;   % Tangent vector in plane.
  v = cross( n, t );

  i = i + 1;
end
v = v / norm( v );
u = cross( v, n );

n_v = size(vertices,1);
vertices_2d = [ dot(vertices,repmat(u,n_v,1),2), ...
                dot(vertices,repmat(v,n_v,1),2) ];

%------------------------------------------------------------------------------%
function [ n, a ] = l_polygon_normal_area( vertices, faces )

if( ~iscell(faces) )
  faces = { faces };
end

n = [];
a = [];
for i=1:length(faces)
  f = faces{i};
  n_e = size(f,2);
  assert( n_e>=3 )

  l_n_i = zeros(size(f,1),n_e);
  a_i = zeros(size(f,1),3);
  for j=1:n_e
    k = mod(j,n_e)+1;
    l = mod(k,n_e)+1;
    p1 = vertices(f(:,j),:);
    p2 = vertices(f(:,k),:);
    p3 = vertices(f(:,l),:);

    n_i = cross( p2-p1, p3-p1, 2 );
    l_n_i(:,j) = sqrt( sum(n_i.^2,2) );
    a_i = a_i + cross( p1, p2, 2 );
  end
  a = [ a; a_i ];

  n_i = zeros(size(f,1),3);
  for jj=1:size(f,1)
    [tmp,j] = max(l_n_i(jj,:));
    k = mod(j,n_e)+1;
    l = mod(k,n_e)+1;
    p1 = vertices(f(jj,j),:);
    p2 = vertices(f(jj,k),:);
    p3 = vertices(f(jj,l),:);

    n_i(jj,:) = cross( p2-p1, p3-p1, 2 );
    n_i(jj,:) = n_i(jj,:)./repmat(l_n_i(jj,j),1,3);
  end
  n = [ n; n_i ];

end
a = dot( n, a, 2 )/2;

%------------------------------------------------------------------------------%
function [ e, eh ] = l_get_edges( faces )

if( ~iscell(faces) )
  faces = { faces };
end

e  = [];   % All edges (n_e x 2).
eh = [];   % Help array where each row consists of face cell group,
% local face, and local edge number (n_e x 3).
for i=1:length(faces)
  f   = faces{i};
  n_f = size(f,1);
  n_e = size(f,2);

  for ie=1:n_e   % Loop over edges.
    iv = ie;   % iv, jv are local vertex numbers for edge ie.
    jv = mod(iv,n_e) + 1;
    e  = [ e; f(:,[iv,jv]) ];

    if( nargout>1 )
      oh = ones(n_f,1);
      eh = [ eh; [i*oh, [1:n_f]', ie*oh] ];
    end
  end
end

%------------------------------------------------------------------------------%
function [ e ] = l_get_external_edges( faces )

e = l_get_edges( faces );
[tmp,ind_e_unq,ixh] = uunique( sort(e,2) );
ind_e_unq = ixh(ind_e_unq);
cnt_e = hist( ind_e_unq, 1:max(ind_e_unq) );
ind_e_ext = find( cnt_e==1 );
e = e(ind_e_ext,:);

%------------------------------------------------------------------------------%
function [ ind_ce ] = l_collinear_edge_vertices( v, e, tol )

if( nargin<3 )
  tol = eps*1e3;
end

n = max(e(:));
ix_ce = zeros(1,n);
for i=1:n
  [i_row,i_col] = find(e==i);
  if( isempty(i_row) )
    continue;
  end
  i_col = 3 - i_col;
  t1  = v(e(i_row(1),i_col(1)),:) - v(i,:);
  t2  = v(e(i_row(2),i_col(2)),:) - v(i,:);
  if( all( abs(cross(t1,t2))<tol ) )
    ix_ce(i) = 1;
  end
end
ind_ce = find(ix_ce);

%------------------------------------------------------------------------------%
function [ vertices, t ] = l_recombine_polygon( polygon, tol )

vertices = deduplicate( polygon{4}, 1, tol );
n_v = size(vertices,1);
if( n_v<=2 )

  vertices = [];
  t = [];

elseif( n_v==3 )

  t  = 1:3;

else

  v2 = l_project_boundary( vertices, {1:n_v} );
  t  = l_delaunay( v2 );

  n1 = l_polygon_normal_area( vertices, 1:n_v );
  n2 = l_polygon_normal_area( vertices, t );
  ix_reorient = dot(repmat(n1,size(n2,1),1),n2,2)<0;
  t(ix_reorient,:) = t(ix_reorient,3:-1:1);

end

%------------------------------------------------------------------------------%
function [ is_in ] = l_in_polygon( p_c, v, f )

n_t = size(p_c,1);
is_in = zeros(n_t,1);
for j=1:length(f)
  f_j = f{j};
  for k=1:size(f_j,1)
    xv = v(f_j(k,:),1);
    yv = v(f_j(k,:),2);
    [in,on] = inpolygon( p_c(:,1), p_c(:,2), xv, yv );
    is_in = is_in | in | on;
  end
end

%------------------------------------------------------------------------------%
function [ t ] = l_delaunay( p )

n_sdim = size(p,2);

if( n_sdim==2 )
  t = delaunay( p(:,1), p(:,2) );
else
  if( exist('DelaunayTri') )
    t = DelaunayTri( p(:,1), p(:,2), p(:,3) );
  elseif( exist('delaunay3') && ~exist( 'OCTAVE_VERSION', 'builtin' ) )
    t = delaunay3( p(:,1), p(:,2), p(:,3) );
  else
    t = delaunay( p );
  end
end
if( isa(t,'DelaunayTri') )
  t = t.Triangulation;
end

%------------------------------------------------------------------------------%
function [ hp ] = plot_polygons( polygons, style, color )
% Plot polygons in CSG tree A.

if( nargin<2 )
  style = 1;
end
if( nargin<3 )
  color = 0.9*[1,1,1];
end
LW = 2;

n_polygons = size(polygons,2);
n_sdim = size(polygons{3,1},2);

hold on
hp = [];
for i=1:n_polygons
  v = polygons{4,i};
  pc = mean(v);

  f = 1:size(v,1);
  hp = [ hp, patch( 'faces', f, 'vertices', v, 'facecolor', color, 'linewidth', LW ) ];


  text( pc(1), pc(2), pc(3), num2str(i) )
  n = l_polygon_normal_area( v, f );
  quiver3(pc(1),pc(2),pc(3),n(1),n(2),n(3))
  if( style==2 )
    set( hp, 'facecolor', 'g' )
    h = plot3( v(:,1), v(:,2), v(:,3), 'r.', 'markersize', 20 );
  end

  if( style==2 )
    % axis equal
    xlabel( 'x ')
    ylabel( 'y ')
    zlabel( 'z ')
    alpha( 0.5 )
    view(3)
    rotate3d('on')
    pause
    set( hp, 'facecolor', color )
    delete( h )
  end
end
% axis equal
xlabel( 'x ')
ylabel( 'y ')
zlabel( 'z ')
try
  alpha( 0.5 )
catch, end
view(3)
rotate3d('on')

if( ~nargout )
  clear hp
end
