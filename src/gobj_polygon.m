function [ gobj ] = gobj_polygon( p, tag )
%GOBJ_POLYGON Create polygon geometry object.
%
%   [ GOBJ ] = GOBJ_POLYGON( P, TAG ) Creates a polygon
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0 0;1 0;0 1]}    Polygon vertex points
%       tag         string {P1}               Geometry object tag/name

% Initial version 150227.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help gobj_polygon, return, end

if( nargin<2 || ~ischar(tag) )
  tag  = 'P1';
end
if( nargin<1 || ~(isnumeric(p) && size(p,1)>=3 && size(p,2)==2) )
  p = [ 0, 0 ;
        1, 0 ;
        0, 1 ];
end


gobj.points = p;
gobj.tag    = tag;
gobj.type   = 'polygon';
gobj.nsdim  = 2;
gobj.v      = p;
gobj.bbox   = [ min(gobj.v); max(gobj.v) ];


% Boundary segments.
n_p = size(p,1);
offset = 1e-2*norm(mean(p));
if( offset<=eps )
  offset = 1e-2*norm(diff(gobj.bbox));
end

for i=1:n_p
  j = mod(i,n_p) + 1;
  boundaries(i).type     = 'line_segment';
  boundaries(i).edges    = [ p(i,:); p(j,:) ];
  boundaries(i).offset   = offset;
  boundaries(i).interior = 0;
  boundaries(i).param    = {};
end
gobj.boundaries = boundaries;


gobj = l_reorient_boundaries( gobj, offset );


%------------------------------------------------------------------------------%
function gobj = l_reorient_boundaries( gobj, offset );
%L_REORIENT_BOUNDARIES Compute vertices and bbox from boundaries.
%
% [ GOBJ ] = L_REORIENT_BOUNDARIES( GOBJ, OFFSET ) Use boundary information
% in GOBJ to reorient boundaries. For edges tests wether the boundary mid point
% shifted by OFFSET (default 1e-6) in the normal direction is inside the geometry.

if( nargin<2 )
  offset = 1e-6;
end

if( isfield(gobj,'boundaries') && ~isempty(gobj.boundaries) )

  boundaries = gobj.boundaries;

  p_inside = zeros(length(boundaries),gobj.nsdim);
  for ib=1:length(boundaries)

    pi   = boundaries(ib).edges;
    i    = ceil(size(pi,1)/2);
    pim  = pi(i:i+1,:);

    pmid = mean(pim);
    t    = pim(2,:) - pim(1,:);
    n    = [ -t(2), t(1) ]/sqrt( t(1)^2 + t(2)^2 );

    if( isfield(boundaries(ib),'offset') )
      p_inside(ib,:) = pmid + boundaries(ib).offset*n;
    else
      p_inside(ib,:) = pmid + offset*n;
    end


  end

  dist = l_compute_distance( gobj, p_inside );

  for ib=1:length(boundaries)
    if( dist(ib)>0 )
      boundaries(ib).edges = boundaries(ib).edges(end:-1:1,:);
    end
  end

  gobj.boundaries = boundaries;

end


%------------------------------------------------------------------------------%
function [ dist ] = l_compute_distance( gobj, p_eval )

n_p = size(p_eval,1);
v   = gobj.points;
n_s = size(v,1);       % Number of segments.
v   = [ v; v(1,:) ];   % Close polygon.

dist = zeros( n_s, n_p );
for i_s=1:n_s
  dist(i_s,:) = l_dist_line( p_eval, v([ i_s i_s+1 ],:) );
end
dist = min( dist )';

% Set negative distances on the inside.
dist = (-1).^( inpolygon( p_eval(:,1), p_eval(:,2), v(:,1), v(:,2) ) ).*dist;

%------------------------------------------------------------------------------%
function [ dist ] = l_dist_line( p, v )

n_p = size(p,1);
w   = v(2,:) - v(1,:);
ix1 = ones(n_p,1);
vp  = v(ix1,:) - p;
w1  = w(ix1,:);
s   = dot( w1, vp, 2);

u = -s/(w*w');
u(u<0) = 0;
u(u>1) = 1;

h = w1.*[u, u] + vp;
dist = sqrt(dot(h,h,2));
