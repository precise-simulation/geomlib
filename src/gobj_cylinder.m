function [ gobj ] = gobj_cylinder( p, r, l, ax, tag, n_s, triangulate )
%GOBJ_CYLINDER Create cylinder geometry object.
%
%   [ GOBJ ] = GOBJ_CYLINDER( P, R, L, AX, TAG, N_S, T ) Creates a cylinder
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0,0,0]}         Coordinates of base center point
%       r           scalar  {1}               Cylinder radius
%       l           scalar  {1}               Cylinder length
%       ax          scalar/array {1}          Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,0])
%       tag         string  {C1}              Geometry object tag/name
%       n_s         scalar  {16}              Number of circumferential boundary segments
%       t           logical {false}           Triangulate boundary segments

% Initial version 150227.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help gobj_cylinder, return, end

if( nargin<7 || ~islogical(triangulate) )
  triangulate = true;
end
if( nargin<6 || ~(isnumeric(n_s) && isscalar(n_s)) )
  n_s = 16;
end
if( nargin<5 || ~ischar(tag) )
  tag = 'C1';
end
if( nargin<4 || ~(isnumeric(ax) && any(length(ax)==[1,3])) )
  ax = 1;
end
if( nargin<3 || ~(isnumeric(l) && isscalar(l)) )
  l = 1;
end
if( nargin<2 || ~(isnumeric(r) && isscalar(r)) )
  r = 1;
end
if( nargin<1 || ~(isnumeric(p) && length(p)==3) )
  p = [0,0,0];
end


gobj.center = p;
gobj.radius = r;
gobj.length = l;
gobj.axis   = ax;
gobj.tag    = tag;
gobj.type   = 'cylinder';
gobj.nsdim  = 3;
gobj.v      = [];
gobj.bbox   = [];

gobj = l_construct_cylinder( gobj, n_s, triangulate );


%------------------------------------------------------------------------------%
function [ gobj ] = l_construct_cylinder( gobj, n_s, triangulate )

assert( rem(n_s/4,1)==0, ...
        'The number of side segments must be subdivisible by 4.' )

c  = gobj.center;
r  = gobj.radius;
l  = gobj.length;
ax = gobj.axis;


% Vertices around z-axis.
n_p = n_s + 1;
th  = linspace( 0, 2*pi, n_p )';
th(end) = th(1);

x = c(1) + r*cos(th);
y = c(2) + r*sin(th);
z = c(3)*ones(n_p,1);
p1 = [ x, y, z ];
p2 = p1;
p2(:,3) = p2(:,3) + l;

% Calculate coordinate rotation.
if( length(ax)==1 )
  switch ax
    case 1
      ax = [1,0,0];
    case 2
      ax = [0,1,0];
    otherwise
      ax = [0,0,1];
  end
elseif( length(ax)~=3 )
  error( 'Not valid axis specification.' )
end
R = l_compute_rotation_matrix( [0,0,1], ax );

% Compute vertices.
p1 = p1 - repmat(c,size(p1,1),1);
p2 = p2 - repmat(c,size(p2,1),1);
p1 = [ R*p1'  ]';
p2 = [ R*p2'  ]';
p1 = p1 + repmat(c,size(p1,1),1);
p2 = p2 + repmat(c,size(p2,1),1);
p  = [ p1; p2 ];

n_ps = n_s/4;
gobj.v = [ p1(1:n_ps:3*n_ps+1,:) ;
           p2(1:n_ps:3*n_ps+1,:) ];
gobj.bbox = [ min(gobj.v); max(gobj.v) ];


% Bottom face.
[boundaries(1:6).type] = deal( 'face' );
if( ~triangulate )
  boundaries(1).faces = [ 1:n_p-1 ];
  boundaries(1).vertices = p1(n_p-1:-1:1,:);
else
  boundaries(1).faces = [ 1:n_p-1 ;
                          2:n_p-1,1 ;
                          n_p*ones(1,n_p-1) ]';
  boundaries(1).vertices = [ p1(n_p-1:-1:1,:) ;
                             c ];
end
boundaries(1).edges = [1:n_p-1;[2:n_p-1,1]]';

% Side faces.
for i=1:4
  ip = [((i-1)*n_ps+1):(i*n_ps+1)];
  f  = [ ip(1:end-1) ;
         ip(2:end)   ;
         ip(2:end)+n_p ;
         ip(1:end-1)+ n_p ]';
  [ix_p,tmp,f] = unique( f );
  f = reshape( f, n_ps, 4 );
  i_bdr = i+1;
  if( ~triangulate )
    boundaries(i_bdr).faces = f;
  else
    boundaries(i_bdr).faces = [ f(:,1:3); f(:,[1,3,4]) ];
  end
  boundaries(i_bdr).vertices = p(ix_p,:);
  boundaries(i_bdr).edges = [ f(:,1), f(:,2) ;
                              f(end,2), f(end,3) ;
                              f(end:-1:1,3), f(end:-1:1,4) ;
                              f(1,4), f(1,1) ];
end

% Top face.
if( ~triangulate )
  boundaries(6).faces = [ 1:n_p-1 ];
  boundaries(6).vertices = p2(1:n_p-1,:);
else
  boundaries(6).faces = [ 1:n_p-1 ;
                          2:n_p-1,1 ;
                          n_p*ones(1,n_p-1) ]';
  boundaries(6).vertices = [ p2(1:n_p-1,:) ;
                             c + [R*[0,0,l]']' ];
end
boundaries(6).edges = [1:n_p-1;[2:n_p-1,1]]';

[boundaries(:).interior] = deal( 0 );

gobj.boundaries = l_reorient_faces( boundaries, mean(gobj.v) );

%------------------------------------------------------------------------------%
function [ boundaries ] = l_reorient_faces( boundaries, p_gobj_center )
% Reorient faces so that normals point outwards.
% Assumes all faces on the boundary are oriented the same way.

for i=1:length(boundaries)
  v = boundaries(i).vertices;
  f = boundaries(i).faces(1,:);

  p_face_center = mean(v(f,:));
  n_face = cross( v(f(2),:)-v(f(1),:), v(f(3),:)-v(f(1),:) );

  if( dot(n_face,p_gobj_center-p_face_center)>0 )
    boundaries(i).faces = boundaries(i).faces(:,end:-1:1);
  end
end

%------------------------------------------------------------------------------%
function [ R ] = l_compute_rotation_matrix( a, b )
% Rotation matrix R for rotating vector a to b.
% ( https://math.stackexchange.com/questions/180418 )

TOL = eps*1e2;

a = a/norm(a);
b = b/norm(b);
if( norm(a-b)<TOL )
  R = eye(3);
  return
end
if( norm(a+b)<TOL )
  R = -eye(3);
  return
end
v = cross(a,b);

ssc = [     0, -v(3),  v(2) ;
         v(3),     0, -v(1) ;
        -v(2),  v(1),     0 ];
R = eye(3) + ssc + ssc^2*(1-dot(a,b))/(norm(v))^2;
