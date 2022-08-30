function [ gobj ] = gobj_circle( p, r, tag )
%GOBJ_CIRCLE Create circle geometry object.
%
%   [ GOBJ ] = GOBJ_CIRCLE( P, R, TAG ) Creates a circle
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0 0]}            Coordinates of center point
%       r           scalar {1}                Circle radius
%       tag         string {C1}               Geometry object tag/name

% Initial version 150227.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help gobj_circle, return, end

if( nargin<3 || ~ischar(tag) )
  tag = 'C1';
end
if( nargin<2 || ~(isnumeric(r) && isscalar(r)) )
  r = 1;
end
if( nargin<1 || ~(isnumeric(p) && length(p)==2) )
  p = [0,0];
end
p = p(:)';


gobj.center = p;
gobj.radius = r;
gobj.tag    = tag;
gobj.type   = 'circle';
gobj.nsdim  = 2;

xcent = p(1);
ycent = p(2);

gobj.v    = [ xcent+r ycent   ;
              xcent   ycent+r ;
              xcent-r ycent   ;
              xcent   ycent-r ];
gobj.bbox = [ min(gobj.v); max(gobj.v) ];


% Approximate boundary with 4*N line segments.
N = 8;
boundary1.type     = 'circle_segment';
boundary1.edges    = circle_line_approx( r, xcent, ycent, 0, pi/2, N );
boundary1.offset   = 1e-2*2*r;
boundary1.interior = 0;
boundary1.param    = { r, xcent, ycent, 0, pi/2 };

boundary2 = boundary1;
boundary2.edges = circle_line_approx( r, xcent, ycent, pi/2, pi, N );
boundary2.param = { r, xcent, ycent, pi/2, pi };

boundary3 = boundary1;
boundary3.edges = circle_line_approx( r, xcent, ycent, pi, 3/2*pi, N );
boundary3.param = { r, xcent, ycent, pi, 3/2*pi };

boundary4 = boundary1;
boundary4.edges = circle_line_approx( r, xcent, ycent, 3/2*pi, 2*pi, N );
boundary4.param = { r, xcent, ycent, 3/2*pi, 2*pi };

gobj.boundaries(1) = boundary1;
gobj.boundaries(2) = boundary2;
gobj.boundaries(3) = boundary3;
gobj.boundaries(4) = boundary4;


%------------------------------------------------------------------------------%
function [ p ] = circle_line_approx( r, xc, yc, th1, th2, N )

th = linspace( th1, th2, N+1 )';
p  = [ r*cos(th)+xc, r*sin(th)+yc ];
