function [ gobj ] = gobj_block( xmin, xmax, ymin, ymax, zmin, zmax, tag, triangulate )
%GOBJ_BLOCK Create block geometry object.
%
%   [ GOBJ ] = GOBJ_BLOCK( XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, TAG, T )
%   Creates a block geometry object. Accepts the following input parameters.
%
%       Parameter   Value/{Default}           Description
%       -----------------------------------------------------------------------------------
%       xmin        scalar  {0}               Minimum x-coordinate
%       xmax        scalar  {1}               Maximum x-coordinate
%       ymin        scalar  {0}               Minimum y-coordinate
%       ymax        scalar  {1}               Maximum y-coordinate
%       zmin        scalar  {0}               Minimum z-coordinate
%       zmax        scalar  {1}               Maximum z-coordinate
%       tag         string  {B1}              Geometry object tag/name
%       t           logical {false}           Triangulate boundary segments

% Initial version 150227.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help gobj_block, return, end

if( nargin<8 || ~islogical(triangulate) )
  triangulate = true;
end
if( nargin<7 || ~ischar(tag) )
  tag = 'B1';
end
if( nargin<6 )
  xmin = 0;
  xmax = 1;
  ymin = 0;
  ymax = 1;
  zmin = 0;
  zmax = 1;
end

tmp  = xmin;
xmin = min( xmin, xmax );
xmax = max( tmp,  xmax );
tmp  = ymin;
ymin = min( ymin, ymax );
ymax = max( tmp,  ymax );
tmp  = zmin;
zmin = min( zmin, zmax );
zmax = max( tmp,  zmax );


gobj.x_min = xmin;
gobj.x_max = xmax;
gobj.y_min = ymin;
gobj.y_max = ymax;
gobj.z_min = zmin;
gobj.z_max = zmax;
gobj.tag   = tag;
gobj.type  = 'block';
gobj.nsdim = 3;


gobj.v = [ xmin, ymin, zmin ;
           xmax, ymin, zmin ;
           xmax, ymax, zmin ;
           xmin, ymax, zmin ;
           xmin, ymin, zmax ;
           xmax, ymin, zmax ;
           xmax, ymax, zmax ;
           xmin, ymax, zmax ];

gobj.bbox = [ min(gobj.v) ;
              max(gobj.v) ];


faces = [ 1, 4, 3, 2 ;
          1, 2, 6, 5 ;
          2, 3, 7, 6 ;
          3, 4, 8, 7 ;
          4, 1, 5, 8 ;
          5, 6, 7, 8 ];

for i_bdr=1:size(faces,1)
  if( ~triangulate )
    faces_i = [ 1:4 ];
  else
    faces_i = [ 1:3; 1,3,4 ];
  end

  boundaries(i_bdr).faces = faces_i;
  boundaries(i_bdr).vertices = gobj.v(faces(i_bdr,:),:);
  boundaries(i_bdr).edges = [ 1:4; 2:4,1 ]';
end

[boundaries(:).interior] = deal( 0 );

gobj.boundaries = boundaries;
