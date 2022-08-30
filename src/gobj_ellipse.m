function [ gobj ] = gobj_ellipse( p, rx, ry, tag )
%GOBJ_ELLIPSE Create ellipse geometry object.
%
%   [ GOBJ ] = GOBJ_ELLIPSE( P, RX, RY, TAG ) Creates an ellipse
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       p           array  {[0 0]}            Coordinates of center point
%       rx          scalar {1}                Radius along x-axis
%       ry          scalar {0.5}              Radius along y-axis
%       tag         string {E1}               Geometry object tag/name

% Initial version 150227.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help gobj_ellipse, return, end

if( nargin<4 || ~ischar(tag) )
  tag  = 'E1';
end
if( nargin<3 || ~(isnumeric(rx) && isscalar(rx)) )
  rx = 1;
end
if( nargin<2 || ~(isnumeric(ry) && isscalar(ry)) )
  ry = 1/2;
end
if( nargin<1 || ~(isnumeric(p) && length(p)==2) )
  p = [0,0];
end
p = p(:)';


gobj.center   = p;
gobj.x_radius = rx;
gobj.y_radius = ry;
gobj.tag      = tag;
gobj.type     = 'ellipse';
gobj.nsdim    = 2;

xcent = p(1);
ycent = p(2);

gobj.v    = [ xcent+rx ycent   ;
              xcent    ycent+ry ;
              xcent-rx ycent   ;
              xcent    ycent-ry ];
gobj.bbox = [ min(gobj.v); max(gobj.v) ];


% Approximate boundary with 4*N line segments.
N = 8;
boundary1.type     = 'ellipse_segment';
boundary1.edges    = ellipse_line_approx( rx, ry, xcent, ycent, 0, pi/2, N );
boundary1.offset   = 1e-2*2*mean([rx ry]);
boundary1.interior = 0;
boundary1.param    = { rx, ry, xcent, ycent, 0, pi/2 };

boundary2 = boundary1;
boundary2.edges = ellipse_line_approx( rx, ry, xcent, ycent, pi/2, pi, N );
boundary2.param = { rx, ry, xcent, ycent, pi/2, pi };

boundary3 = boundary1;
boundary3.edges = ellipse_line_approx( rx, ry, xcent, ycent, pi, 3/2*pi, N );
boundary3.param = { rx, ry, xcent, ycent, pi, 3/2*pi };

boundary4 = boundary1;
boundary4.edges = ellipse_line_approx( rx, ry, xcent, ycent, 3/2*pi, 2*pi, N );
boundary4.param = { rx, ry, xcent, ycent, 3/2*pi, 2*pi };

gobj.boundaries(1) = boundary1;
gobj.boundaries(2) = boundary2;
gobj.boundaries(3) = boundary3;
gobj.boundaries(4) = boundary4;


%------------------------------------------------------------------------------%
function [ p ] = ellipse_line_approx( rx, ry, xc, yc, th1, th2, N )

  th = linspace( th1, th2, N+1 )';
  p  = [ rx*cos(th)+xc, ry*sin(th)+yc ];
