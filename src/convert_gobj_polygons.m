function [ polygons ] = convert_gobj_polygons( gobj, id )
%CONVERT_GOBJ_POLYGONS Constructs polygons from geometry object.
%
%   [ POLYGONS ] = CONVERT_GOBJ_POLYGONS( GOBJ, ID ) Constructs
%   polygons from the boundaries of geometry object GOBJ by calling
%   the CSG_OP build operation. ID if present will be appended to the
%   polygon id property. Thus the final polygon identity field will
%   consist of the local boundary number, local polygon number, and ID.

% Initial version 171228.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.
if( ~(nargin || nargout) ),help convert_gobj_polygons, return, end

if( nargin<2 || ~isnumeric(id) )
  id = 0;
end


if( isfield(gobj,'boundaries') )

  polygons = {};
  boundaries = gobj.boundaries;
  for i_bdr=1:length(boundaries)
    polygons_ibdr = boundary_to_polygons( boundaries(i_bdr), [i_bdr,id] );
    polygons = [ polygons, polygons_ibdr ];
  end

elseif( isfield(gobj,'edges') || isfield(gobj,'faces') )   % Input is boundary struct.

  polygons = boundary_to_polygons( gobj, id );

else
  error( 'Could not create polygons from geometry object.' )
end

%------------------------------------------------------------------------------%
function [ polygons ] = boundary_to_polygons( boundary, id )

polygons = {};
if( isfield(boundary,'vertices') && isfield(boundary,'faces') )

  vertices = boundary.vertices;
  faces = boundary.faces;
  if( ~iscell(faces) )
    faces = { faces };
  end

  i_cnt = 0;
  for i=1:length(faces)
    f = faces{i};
    for j=1:size(f,1)
      v = vertices(f(j,:),:);
      i_cnt = i_cnt + 1;
      identity = [i_cnt,id];

      polygon_i = csg_op( v, identity, 'b' );
      polygons  = [ polygons, polygon_i ];
    end
  end

elseif( isfield(boundary,'edges') )

  p_e = boundary.edges;
  n_e = size(p_e,1) - 1;
  for i=1:n_e
    v = p_e([i,i+1],:);
    identity = [i,id];

    polygon_i = csg_op( v, identity, 'b' );
    polygons  = [ polygons, polygon_i ];
  end

else
  error( 'Could not create polygons from boundary.' )
end
