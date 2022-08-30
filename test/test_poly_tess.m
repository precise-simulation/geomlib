function [ res ] = test_poly_tess( fun )
if(nargin),if(ischar(fun)),fun=str2func(fun);end,fun();if(nargout>0),res=[];end,return,end
res = mloct_test_caller(cellfun(@str2func,get_subfuns([mfilename('fullpath'),'.m'],'test'),'UniformOutput',0));


function test_empty

[vertices,faces,ids] = csg_polygon_tesselation( [], 1 );

assert( isempty(vertices) )
assert( isempty(faces) )
assert( isempty(ids) )


function test_tri

v = [0 0 0;1 0 0;0 1 0];
polygon = csg_op( v, 1, 'b' );
[vertices,faces,ids] = csg_polygon_tesselation( polygon );

assert( isequal(v,vertices) )
assert( isequal(1:3,faces) )
assert( ids==1 )


function test_quad

v = [0 0 0;1 0 0;1 1 0;0 1 0];
polygon = csg_op( v, 2, 'b' );
[vertices,faces,ids] = csg_polygon_tesselation( polygon );

assert( isequal(v,vertices) )
assert( isequal([3 4 2;4 1 2],faces) || ...
        isequal([1 2 4;4 2 3],faces) || ...
        isequal([1 2 3;4 1 3],faces) )
assert( isequal(ids,[2;2]) )


function test_lshape

v1 = [0 0 0;1 0 0;1 1 0;0 2 0];
v2 = [1 1 0;2 1 0;2 2 0;0 2 0];
polygons = [ csg_op( v1, [3 4], 'b' ), ...
             csg_op( v2, [3 4], 'b' ) ];
[vertices,faces,ids] = csg_polygon_tesselation( polygons );

assert( size(vertices,1)==6 )
assert( size(faces,1)==4 )
assert( isequal(ids,[3 4;3 4;3 4;3 4]) )
for i=1:4
  assert( inpolygon(1.5,0.5,vertices(faces(i,:),1),vertices(faces(i,:),2))==0 )
end


function test_no_merge

v1 = [0 0 0;1 0 0;1 1 0;0 2 0];
v2 = [1 1 0;2 1 0;2 2 0;0 2 0];
polygons = [ csg_op( v1, 5, 'b' ), ...
             csg_op( v2, 6, 'b' ) ];
[vertices,faces,ids] = csg_polygon_tesselation( polygons );

assert( isequal([v1;v2],vertices ) )
assert( size(faces,1)==4 )
assert( isequal(ids,[5;5;6;6]) )
for i=1:4
  assert( inpolygon(1.5,0.5,vertices(faces(i,:),1),vertices(faces(i,:),2))==0 )
end


function test_patch

load csg_node

[vertices,faces,ids] = csg_polygon_tesselation( csg_node{5} );

assert( size(vertices,1)==4 )
      assert( size(faces,1)==2 )
      assert( isequal(ids,[6 4 2;6 4 2]) )
