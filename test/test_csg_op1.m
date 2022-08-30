function [ res ] = test_csg_op1( fun )
if(nargin),if(ischar(fun)),fun=str2func(fun);end,fun();if(nargout>0),res=[];end,return,end
res = mloct_test_caller(cellfun(@str2func,get_subfuns([mfilename('fullpath'),'.m'],'test'),'UniformOutput',0));


function test_non_intersecting_rectangles

for s = [ 1e-6, 1e-3, 1, 1e3 ]
  xmin = s*1;
  xmax = s*2;
  ymin = s*3;
  ymax = s*4;

  o1 = gobj_rectangle(xmin,xmax,ymin,ymax);
  p1 = convert_gobj_polygons( o1 );
  a  = csg_op( p1, 'b' );

  for dh = [ 1e-6, 1e-3, 1 ]
    for i=1:3
      x1 = s*(i-1) + (i-2)*dh;
      x2 = s*i     + (i-2)*dh;
      for j=1:3
        y1 = s*(j+1) + (j-2)*dh;
        y2 = s*(j+2) + (j-2)*dh;
        if(i==j),continue,end

        o2 = gobj_rectangle(x1,x2,y1,y2);
        p2 = convert_gobj_polygons( o2 );
        b  = csg_op( p2, 'b' );

        s_msg = sprintf( 'Assertion 1 failed: s = %g, dh = %g, ij = [%i,%i]', s, dh, i, j );
        [c,d,stat] = csg_op( a, b, '-' );
        assert( stat>0 && isequal(a,c) && isequal(b,d), s_msg )

        s_msg = sprintf( 'Assertion 2 failed: s = %g, dh = %g, ij = [%i,%i]', s, dh, i, j );
        [c,d,stat] = csg_op( b, a, '-' );
        assert( stat>0 && isequal(a,d) && isequal(b,c), s_msg )

        s_msg = sprintf( 'Assertion 3 failed: s = %g, dh = %g, ij = [%i,%i]', s, dh, i, j );
        [c,d,stat] = csg_op( a, b, '&' );
        assert( stat>0 && isequal(a,c) && isequal(b,d), s_msg )

      end
    end
  end

end


function test_non_intersecting_blocks

for s = [ 1e-6, 1e-3, 1, 1e3 ]
  xmin = s*1;
  xmax = s*2;
  ymin = s*3;
  ymax = s*4;
  zmin = s*5;
  zmax = s*6;

  o1 = gobj_block(xmin,xmax,ymin,ymax,zmin,zmax);
  p1 = convert_gobj_polygons( o1 );
  a  = csg_op( p1, 'b' );

  for dh = [ 1e-6, 1e-3, 1 ]
    for i=1:3
      x1 = s*(i-1) + (i-2)*dh;
      x2 = s*i     + (i-2)*dh;
      for j=1:3
        y1 = s*(j+1) + (j-2)*dh;
        y2 = s*(j+2) + (j-2)*dh;
        for k=1:3
          z1 = s*(k+3) + (k-2)*dh;
          z2 = s*(k+4) + (k-2)*dh;
          if(i==j && j==k),continue,end

          o2 = gobj_block(x1,x2,y1,y2,z1,z2);
          p2 = convert_gobj_polygons( o2 );
          b  = csg_op( p2, 'b' );

          s_msg = sprintf( 'Assertion 1 failed: s = %g, dh = %g, ij = [%i,%i]', s, dh, i, j );
          [c,d,stat] = csg_op( a, b, '-' );
          assert( stat>0 && isequal(a,c) && isequal(b,d), s_msg )

          s_msg = sprintf( 'Assertion 2 failed: s = %g, dh = %g, ij = [%i,%i]', s, dh, i, j );
          [c,d,stat] = csg_op( b, a, '-' );
          assert( stat>0 && isequal(a,d) && isequal(b,c), s_msg )

          s_msg = sprintf( 'Assertion 3 failed: s = %g, dh = %g, ij = [%i,%i]', s, dh, i, j );
          [c,d,stat] = csg_op( a, b, '&' );
          assert( stat>0 && isequal(a,c) && isequal(b,d), s_msg )

        end
      end
    end
  end

end


function test_intersecting_rectangles

for s = [ 1e-5, 1e-3, 1, 1e3 ]
  xmin = s*1;
  xmax = s*2;
  ymin = s*3;
  ymax = s*4;

  o1 = gobj_rectangle(xmin,xmax,ymin,ymax);
  p1 = convert_gobj_polygons( o1 );
  a  = csg_op( p1, 'b' );

  for dx = [ -s/10, -s/2, s/2, s/10 ]
    for dy = [ -s/10, -s/2, s/2, s/10 ]

      o2 = gobj_rectangle(xmin+dx,xmax+dx,ymin+dy,ymax+dy);
      p2 = convert_gobj_polygons( o2 );
      b  = csg_op( p2, 'b' );


      s_msg = sprintf( 'Assertion 1 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
      [c,d,stat] = csg_op( a, b, '+' );
      p = csg_op( c, 'e' );
      n_p = size(p,2);
      assert( stat==0 && n_p>4 && n_p<=10 && isempty(d), s_msg )

      s_msg = sprintf( 'Assertion 2 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
      [c,d,stat] = csg_op( a, b, '-' );
      p = csg_op( c, 'e' );
      n_p = size(p,2);
      assert( stat==0 && n_p>4 && n_p<=8 && isempty(d), s_msg )

      s_msg = sprintf( 'Assertion 3 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
      [c,d,stat] = csg_op( b, a, '-' );
      p = csg_op( c, 'e' );
      n_p = size(p,2);
      assert( stat==0 && n_p>4 && n_p<=8 && isempty(d), s_msg )

      s_msg = sprintf( 'Assertion 4 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
      [c,d,stat] = csg_op( b, a, '&' );
      p = csg_op( c, 'e' );
      n_p = size(p,2);
      assert( stat==0 && n_p==4 && isempty(d), s_msg )

    end
  end

end


function test_intersecting_blocks

for s = [ 1e-3, 1, 1e3 ]
  xmin = s*1;
  xmax = s*2;
  ymin = s*3;
  ymax = s*4;
  zmin = s*5;
  zmax = s*6;

  o1 = gobj_block(xmin,xmax,ymin,ymax,zmin,zmax);
  p1 = convert_gobj_polygons( o1 );
  a  = csg_op( p1, 'b' );

  for dx = [ -s/10, -s/2, s/2, s/10 ]
    for dy = [ -s/10, -s/2, s/2, s/10 ]
      for dz = [ -s/10, -s/2, s/2, s/10 ]

        o2 = gobj_block(xmin+dx,xmax+dx,ymin+dy,ymax+dy,zmin+dz,zmax+dz);
        p2 = convert_gobj_polygons( o2 );
        b  = csg_op( p2, 'b' );

        s_msg = sprintf( 'Assertion 1 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
        [c,d,stat] = csg_op( a, b, '+' );
        p = csg_op( c, 'e' );
        n_p = size(p,2);
        assert( stat==0 && n_p>6 && n_p<=43 && isempty(d), s_msg )

        s_msg = sprintf( 'Assertion 2 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
        [c,d,stat] = csg_op( a, b, '-' );
        p = csg_op( c, 'e' );
        n_p = size(p,2);
        assert( stat==0 && n_p>6 && n_p<=43 && isempty(d), s_msg )

        s_msg = sprintf( 'Assertion 3 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
        [c,d,stat] = csg_op( b, a, '-' );
        p = csg_op( c, 'e' );
        n_p = size(p,2);
        assert( stat==0 && n_p>6 && n_p<=43 && isempty(d), s_msg )

        s_msg = sprintf( 'Assertion 4 failed: s = %g, dx = %g, dy = %g', s, dx, dy );
        [c,d,stat] = csg_op( b, a, '&' );
        p = csg_op( c, 'e' );
        n_p = size(p,2);
        assert( stat==0 && n_p>=8 && n_p<=12 && isempty(d), s_msg )

      end
    end
  end

end
