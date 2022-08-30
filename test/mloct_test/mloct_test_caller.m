function [ argout ] = mloct_test_caller( subFcns )
%MLOCT_TEST_CALLER Test caller function
%
% [ ARGOUT ] = MLOCT_TEST_CALLER( SUBFCNS ) Test caller driver function.
% Called from the test suite function with test function handles SUBFCNS
% as input. Returns a cell array with the NPASSFAIL array and cell array
% of error message structs.

% Initial version 180216.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.

[fid,isVerbose,nCnt] = mloct_test_setup();
nSplit = 20;

% Setup.
[ST,I] = dbstack( '-completenames' );
sCallerName = ST(I+1).name;
sCallerFile = ST(I+1).file;
if( nargin<1 )
  subFcns = which( '-subfun', sCallerFile );
end
if( isempty(subFcns) )
  subFcns = get_subfuns( sCallerFile, 'test' );
end


sSuite = sCallerName;
cTests = subFcns;
for i=length(subFcns):-1:1
  sFcn = subFcns{i};
  if( isa(sFcn,'function_handle') )
    sFcn = functions(sFcn);
    sFcn = sFcn.function;
  end
  if( ~strncmpi(sFcn,'test',4) )
    subFcns(i) = [];
    cTests(i)  = [];
  else
    cTests{i} = sFcn;
  end
end
if( ~strncmpi(sCallerName,'test',4) || isempty(cTests) )
  return
end


% Test loop.
nPassFail = [0,0];   % Number of passed/failed tests.
cErrors = {};
tTot = 0;
for i=1:length(cTests)
  sTest = cTests{i};

  if( isVerbose )
    fprintf(fid,'   %s %s ',sTest,repmat('.',1,57-length(sTest)));
  end

  tic
  try

    if( isa(subFcns{i},'function_handle') )
      subFcns{i}();
    else
      evalin( 'caller', sTest )
    end
    t = toc;

    if( isVerbose )
      fprintf(fid,'passed in %12.6f seconds\n', t);
    else
      fprintf(fid,'.', t);
      if( mod(nCnt+i,nSplit)==0 )
        fprintf(fid,'\n');
      end
    end
    nPassFail(1) = nPassFail(1) + 1;

  catch me

    t = toc;
    msg.message = me.message;
    if( isfield(me,'stack') )
      msg.stack = me.stack;
    else
      msg.stack.file = sCallerName;
      msg.stack.name = sTest;
      msg.stack.line = [];
    end
    msg.sCallerName = sCallerName;
    msg.sTest = sTest;
    cErrors = [ cErrors, {msg} ];

    if( isVerbose )
      fprintf(fid,'FAILED in %12.6f seconds\n', t);
    else
      fprintf(fid,'F', t);
      if( mod(nCnt+i,nSplit)==0 )
        fprintf(fid,'\n');
      end
    end
    nPassFail(2) = nPassFail(2) + 1;

  end
  drawnow
  tTot = tTot + t;

end
if( isVerbose )
  fprintf(fid,'%s ',sSuite,repmat('.',1,60-length(sSuite)));
  if( nPassFail(2)==0 )
    fprintf(fid,'passed in %12.6f seconds\n', tTot);
  else
    fprintf(fid,'FAILED in %12.6f seconds\n', tTot);
  end
end
mloct_test_setup( [], [], nCnt+length(cTests) );

if( nargout )
  argout = { nPassFail, cErrors };
end
