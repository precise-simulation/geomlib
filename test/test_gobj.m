function [ res ] = test_gobj( fun )
if(nargin),if(ischar(fun)),fun=str2func(fun);end,fun();if(nargout>0),res=[];end,return,end
res = mloct_test_caller(cellfun(@str2func,get_subfuns([mfilename('fullpath'),'.m'],'test'),'UniformOutput',0));


function test_gobj_rectangle

o = gobj_rectangle(1,2,3,4);

or.x_min = 1;
or.x_max = 2;
or.y_min = 3;
or.y_max = 4;
or.tag   = 'R1';
or.type  = 'rectangle';
or.nsdim = 2;
or.v = [ 1, 3; 2, 3; 2, 4; 1, 4 ];
or.bbox = [min(or.v);max(or.v)];

o = rmfield(o,'boundaries');
assert( isequal(o,or) )


function test_gobj_polygon

o = gobj_polygon([0,0;1,0;0,1]);

or.points = [ 0, 0; 1, 0; 0, 1 ];
or.tag = 'P1';
or.type = 'polygon';
or.nsdim = 2;
or.v = [ 0, 0; 1, 0; 0, 1 ];
or.bbox = [min(or.v);max(or.v)];

o = rmfield(o,'boundaries');
assert( isequal(o,or) )


function test_gobj_circle

p = rand(1,2);
r = rand(1);
o = gobj_circle(p,r);

or.center = p;
or.radius = r;
or.tag = 'C1';
or.type = 'circle';
or.nsdim = 2;
or.v = [ p(1)+r, p(2)   ;
         p(1),   p(2)+r ;
         p(1)-r, p(2)   ;
         p(1),   p(2)-r ];
or.bbox = [min(or.v);max(or.v)];

o = rmfield(o,'boundaries');
assert( isequal(o,or) )


function test_gobj_ellipse

p = rand(1,2);
r = rand(1,2);
o = gobj_ellipse(p,r(1),r(2));

or.center = p;
or.x_radius = r(1);
or.y_radius = r(2);
or.tag = 'E1';
or.type = 'ellipse';
or.nsdim = 2;
or.v = [ p(1)+r(1), p(2)   ;
         p(1),      p(2)+r(2) ;
         p(1)-r(1), p(2)   ;
         p(1),      p(2)-r(2) ];
or.bbox = [min(or.v);max(or.v)];

o = rmfield(o,'boundaries');
assert( isequal(o,or) )


function test_gobj_block

o = gobj_block(1,2,3,4,5,6);

or.x_min = 1;
or.x_max = 2;
or.y_min = 3;
or.y_max = 4;
or.z_min = 5;
or.z_max = 6;
or.tag   = 'B1';
or.type  = 'block';
or.nsdim = 3;
or.v = [ 1, 3, 5; 2, 3, 5; 2, 4, 5; 1, 4, 5;
         1, 3, 6; 2, 3, 6; 2, 4, 6; 1, 4, 6 ];
or.bbox = [min(or.v);max(or.v)];

o = rmfield(o,'boundaries');
assert( isequal(o,or) )


function test_gobj_cylinder

p = rand(1,3);
r = rand(1);
l = rand(1);
o = gobj_cylinder(p,r,l,3);

or.center = p;
or.radius = r;
or.length = l;
or.axis = 3;
or.tag = 'C1';
or.type = 'cylinder';
or.nsdim = 3;
v = [ p(1)+r, p(2)  , p(3) ;
      p(1),   p(2)+r, p(3) ;
      p(1)-r, p(2)  , p(3) ;
      p(1),   p(2)-r, p(3) ;
      p(1)+r, p(2)  , p(3)+l ;
      p(1),   p(2)+r, p(3)+l ;
      p(1)-r, p(2)  , p(3)+l ;
      p(1),   p(2)-r, p(3)+l ];

assert( norm(o.v(:)-v(:))<eps*1e3 )

o = rmfield(o,'v');
o = rmfield(o,'bbox');
o = rmfield(o,'boundaries');
assert( isequal(o,or) )


function test_gobj_sphere

p = rand(1,3);
r = rand(1);
o = gobj_sphere(p,r,3);

or.center = p;
or.radius = r;
or.axis = 3;
or.tag = 'S1';
or.type = 'sphere';
or.nsdim = 3;
or.v = [ p(1),   p(2),   p(3)-r ;
         p(1),   p(2)-r, p(3) ;
         p(1)+r, p(2),   p(3) ;
         p(1),   p(2)+r, p(3) ;
         p(1)-r, p(2),   p(3) ;
         p(1),   p(2),   p(3)+r];
or.bbox = [min(or.v);max(or.v)];

o = rmfield(o,'boundaries');
assert( isequal(o,or) )
