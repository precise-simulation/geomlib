function [ gobj ] = gobj_rectangle( xmin, xmax, ymin, ymax, tag )
%GOBJ_RECTANGLE Create rectangle geometry object.
%
%   [ GOBJ ] = GOBJ_RECTANGLE( XMIN, XMAX, YMIN, YMAX, TAG ) Creates a rectangle
%   geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       xmin        scalar {0}                Minimum x-coordinate
%       xmax        scalar {1}                Maximum x-coordinate
%       ymin        scalar {0}                Minimum y-coordinate
%       ymax        scalar {1}                Maximum y-coordinate
%       tag         string {R1}               Geometry object tag/name

% Initial version 150227.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help gobj_rectangle, return, end

if( nargin<5 || ~ischar(tag) )
  tag  = 'R1';
end
if( nargin<4 )
  xmin = 0;
  xmax = 1;
  ymin = 0;
  ymax = 1;
end

tmp  = xmin;
xmin = min( xmin, xmax );
xmax = max( tmp,  xmax );
tmp  = ymin;
ymin = min( ymin, ymax );
ymax = max( tmp,  ymax );


gobj.x_min = xmin;
gobj.x_max = xmax;
gobj.y_min = ymin;
gobj.y_max = ymax;
gobj.tag   = tag;
gobj.type  = 'rectangle';
gobj.nsdim = 2;

gobj.v = [ xmin ymin ;
           xmax ymin ;
           xmax ymax ;
           xmin ymax ];

gobj.bbox = [ min(gobj.v) ;
              max(gobj.v) ];


% Boundary segments.
boundary1.type     = 'line_segment';
boundary1.edges    = [xmin ymin; xmax ymin];
boundary1.offset   = 1e-2*mean( [ xmax-xmin ymax-ymin ] );
boundary1.interior = 0;
boundary1.param    = {};

boundary2 = boundary1;
boundary2.edges = [xmax ymin; xmax ymax];

boundary3 = boundary1;
boundary3.edges = [xmax ymax; xmin ymax];

boundary4 = boundary1;
boundary4.edges = [xmin ymax; xmin ymin];

gobj.boundaries(1) = boundary1;
gobj.boundaries(2) = boundary2;
gobj.boundaries(3) = boundary3;
gobj.boundaries(4) = boundary4;
