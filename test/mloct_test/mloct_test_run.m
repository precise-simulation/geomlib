function [ results ] = mloct_test_run( sSrcDir, sTestDir, sTest, fid, isVerbose )
%MLOCT_TEST_RUN Test run function
%
% [ RESULTS ] = MLOCT_TEST_RUN( SSRCDIR, STESTDIR, STEST, FID, ISVERBOSE )
% Runs test cases in STESTDIR and adds SSRCDIR to the path. STEST can be
% a specific test case function, and FID a file handle. ISVERBOSE toggles
% verbose or short output format. RESULTS is a NTESTSUITES x 2 array with
% the number of passed tests in first column and failed in second.

% Initial version 180216.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.


t1 = tic;
if( nargin )
  sSrcDir  = 'src';
end
addpath( genpath(sSrcDir) );
addpath( genpath(sTestDir) );
mloct_test_setup( fid, isVerbose, 0 );


% Find test (suite) files to run.
sCase = '';
indCase = find(sTest== ':',1);
if( ~isempty(indCase) )
  sCase = sTest(indCase+1:end);
  sTest = sTest(1:indCase-1);
end
if( ~strncmp('test',sTest,4) )
  sTest = ['test',sTest];
end
cTestFiles = struct2cell( dir(sTestDir) );
ind = strmatch( sTest, cTestFiles(1,:) );
if( isempty(ind) )
  sTest = ['test_',sTest(5:end)];
  ind = strmatch( sTest, cTestFiles(1,:) );
end
cTestFiles = cTestFiles(1,ind);
if( isempty(cTestFiles) )
  warning( 'No tests to run.' )
  results = [];
  return
end


% Header.
sMatlab = 'Matlab';
if( exist('OCTAVE_VERSION','builtin') )
  sMatlab = 'Octave';
end
fprintf(fid,'\nRunning test suites on %s %s %s:\n\n', ...
        sMatlab, version(), datestr(clock) );


% Run tests.
results = {};
for i=1:length(cTestFiles)
  sTestFcn = cTestFiles{i};
  if( length(sTestFcn)>2 && strmatch(sTestFcn(end-1:end),'.m') )
    sTestFcn = sTestFcn(1:end-2);
    if( isVerbose )
      fprintf(fid,'%s\n',sTestFcn);
    end

    if( isempty(sCase) )
      res_i = feval( sTestFcn );
    else
      fprintf(fid,'   %s %s ',sCase,repmat('.',1,57-length(sCase)));
      tic
      try
        feval( sTestFcn, sCase );
        res_i = {[1,0],[]};
        fprintf(fid,'passed in %12.6f seconds\n', toc );
      catch me
        fprintf(fid,'FAILED in %12.6f seconds\n', toc );
        msg.message = me.message;
        if( isfield(me,'stack') )
          msg.stack = me.stack;
        else
          msg.stack.file = sTestFcn;
          msg.stack.name = '';
          msg.stack.line = [];
        end
        msg.sCallerName = sTestFcn;
        msg.sTest = sCase;
        res_i = {[0,1],{msg}};
      end

    end
    results = [ results; res_i ];

    if( isVerbose )
      fprintf( fid, '\n' );
    end
  end
end
if( ~isVerbose )
  fprintf( fid, '\n\n' );
end

tTot = toc(t1);

cErrors = results(:,2);
results = vertcat( results{:,1} );

if( ~isempty(cErrors) )
  l_proc_errors( fid, cErrors )
  fprintf( fid, '\n' );
end

fprintf( fid, '\n' );
if( ~any(results(:,2)) )
  fprintf(fid,'Finished test runs with all (%i) PASSED tests in %f seconds.\n\n', sum(results(:)), tTot );
else
  fprintf(fid,'Finished test runs with %i (%i) FAILED tests in %f seconds.\n\n', ...
          sum(results(:,2)), sum(results(:)), tTot );
end

if( ~nargout )
  clear results
end


%------------------------------------------------------------------------------%
function l_proc_errors( fid, cErrors, type )

if( nargin<3 )
  type = 'extended';
end
dash = '-';

for i=1:length(cErrors)
  for j=1:length(cErrors{i})
    err = cErrors{i}{j};

    if( strcmpi(type,'basic') )

      report = [err.message];

    else

      report = sprintf( '%s\nTest case [ %s:%s ] failed\n%s\n', ...
                        dash(ones(1,80)), err.sCallerName, err.sTest, dash(ones(1,80)) );

      for k=1:length(err.stack)
        [tmp1,file,tmp2] = fileparts(err.stack(k).file);
        name = err.stack(k).name;
        line = err.stack(k).line;

        if( strcmpi(name,'mloct_test_caller') )
          break
        end

        s_line = ['Error in ',file,'>',name,' (line ',num2str(line),')'];
        report = sprintf( '  %s\n', report, s_line );
        if( k==1 )
          report = sprintf( '%s    %s\n', report, err.message );
        end

      end
    end

    report = strtrim(report);
    fprintf( fid, '\n%s\n', report );
  end
end
