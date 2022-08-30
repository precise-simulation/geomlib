function [ gobj ] = gobj_sphere( p, r, ax, tag, n_s, triangulate )
%GOBJ_SPHERE Create sphere geometry object.
%
%   [ GOBJ ] = GOBJ_SPHERE( P, R, AX, TAG, N_S, T ) Creates a sphere
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array   {[0 0 0]}         Coordinates of center point
%       r           scalar  {1}               Sphere radius
%       ax          scalar/array {1}          Axis direction (1/2/3 = x/y/z-axis)
%                                             alt. axis direction vector (ex. [1,1,0])
%       tag         string  {S1}              Geometry object tag/name
%       n_s         scalar  {16}              Number of circumferential boundary segments
%       t           logical {false}           Triangulate boundary segments

% Initial version 150227.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help gobj_sphere, return, end

if( nargin<6 || ~islogical(triangulate) )
  triangulate = true;
end
if( nargin<5 || ~(isnumeric(n_s) && isscalar(n_s)) )
  n_s = 16;
end
if( nargin<4 || ~ischar(tag) )
  tag = 'S1';
end
if( nargin<3 || ~(isnumeric(ax) && any(length(ax)==[1,3])) )
  ax = 1;
end
if( nargin<2 || ~(isnumeric(r) && isscalar(r)) )
  r = 1;
end
if( nargin<1 || ~(isnumeric(p) && length(p)==3) )
  p = [0,0,0];
end


gobj.center = p;
gobj.radius = r;
gobj.axis   = ax;
gobj.tag    = tag;
gobj.type   = 'sphere';
gobj.nsdim  = 3;
gobj.v      = [];
gobj.bbox   = [];

gobj = l_construct_sphere( gobj, n_s, triangulate );


%------------------------------------------------------------------------------%
function [ gobj ] = l_construct_sphere( gobj, n_s, triangulate )

assert( rem(n_s/4,1)==0, ...
        'The number of side segments must be subdivisible by 4.' )

c  = gobj.center;
r  = gobj.radius;
ax = gobj.axis;

% Calculate coordinate rotation matrix.
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

% Vertices around z-axis.
n_p = n_s + 1;
for i=1:2     % Bottom/top half (ph).
  if( i==1 )
    ph = linspace( -pi/2, 0, n_s/2+1 );
  else
    ph = linspace( 0, pi/2, n_s/2+1 );
  end

  for j=1:4   % x/y-quadrant (th).
    th = linspace( (j-1)*pi/2, j*pi/2, n_s/4+1 );   % Angles in x/y-plane.

    [ph_ij,th_ij] = meshgrid(ph,th);

    x = c(1) + r*cos(ph_ij).*cos(th_ij);
    y = c(2) + r*cos(ph_ij).*sin(th_ij);
    z = c(3) + r*sin(ph_ij);
    p = [x(:),y(:),z(:)];
    p = p - repmat(c,size(p,1),1);
    p = [ R*p' ]';
    p = p + repmat(c,size(p,1),1);

    ind = reshape(1:numel(x),size(x));
    ind = reshape(ind(1:n_s/4,1:n_s/2),1,n_s^2/8);
    f = [ ind             ;
          ind         + 1 ;
          ind + n_s/4 + 2 ;
          ind + n_s/4 + 1 ]';

    % Remove multiple bottom/top points.
    if( i==1 )
      p(2:n_s/4+1,:) = [];
      ix = ismember(f(:)',2:(n_s/4+1));
      f(ix) = 1;
      ix = f>1;
      f(ix) = f(ix) - n_s/4;

      fti = f(1:n_s/4,2:4);
      fq  = f(n_s/4+1:end,:);
    else
      n_p = (n_s/4+1)*(n_s/2+1);
      p(n_p-n_s/4+1:n_p,:) = [];
      ix = ismember(f(:)',(n_p-n_s/4+1):n_p);
      f(ix) = n_p - (n_s/4);

      fti = f(end-n_s/4+1:end,1:3);
      fq  = f(1:end-n_s/4,:);
    end

    ft = [ fti; fq(:,1:3); fq(:,[1,3,4]) ];
    if( ~triangulate )
      f = { fti; fq };
    else
      f = { ft };
    end

    i_bdr = 4*(i-1) + j;
    boundaries(i_bdr).faces = f;
    boundaries(i_bdr).vertices = p;
    boundaries(i_bdr).edges = l_construct_edges( ft );
  end
end

p_v = [0,0,-r; 0,-r,0; r,0,0; 0,r,0; -r,0,0; 0,0,r];
p_v = repmat( c, 6, 1 ) + [ R*p_v' ]';
gobj.v = p_v;
gobj.bbox = [ min(gobj.v); max(gobj.v) ];

[boundaries(:).interior] = deal( 0 );

gobj.boundaries = l_reorient_faces( boundaries, mean(gobj.v) );


if( triangulate )   % Convert cell to array.
  for i_bdr=1:length(boundaries)
    faces = vertcat( gobj.boundaries(i_bdr).faces{:} );
    gobj.boundaries(i_bdr).faces = faces;
  end
end

%------------------------------------------------------------------------------%
function [ boundaries ] = l_reorient_faces( boundaries, p_gobj_center )
% Reorient faces so that normals point outwards.
% Assumes all faces on the boundary are oriented the same way.

for i=1:length(boundaries)
  v = boundaries(i).vertices;

  f_i = boundaries(i).faces;
  for j=1:length(f_i)
    f = f_i{j}(1,:);

    p_face_center = mean(v(f,:));
    n_face = cross( v(f(2),:)-v(f(1),:), v(f(3),:)-v(f(1),:) );

    if( dot(n_face,p_gobj_center-p_face_center)>0 )
      boundaries(i).faces{j} = boundaries(i).faces{j}(:,end:-1:1);
    end
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

%------------------------------------------------------------------------------%
function [ e ] = l_construct_edges( f )
% Constructs boundary edges from faces.

assert( isnumeric(f) & size(f,2)==3, 'All faces must be triangles.' )

e = [ f(:,1:2); f(:,2:3); f(:,[3,1]) ];
e = sort(e,2);
[e,~,ix] = unique(e,'rows');
cnt = hist(ix,1:max(ix));
e = e( find(cnt==1), : );
