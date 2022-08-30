function [ results ] = runtests( varargin )
%RUNTESTS Runs the testsuite
%
% [ RESULTS ] = RUNTESTS( VARARGIN ) Runs the test suite (VERBOSELY if
% any input argument is given). To run a single test suite/case enter
% the testSuite:testCase can be specified as an argument. If two input
% arguments are given verbose tests are run and output to the file
% given by the second argument. RESULTS is a NTESTSUITES x 2 array
% with the number of passed tests in first column and failed in second.

% The SSRCDIR with source files (which are added to the path) and
% STESTDIR with mloct_test files are specified manually below.

% Initial version 180216.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.

sSrcDir   = 'src';
sTestDir  = 'test';
sTest     = '';
fid       = 1;
isVerbose = false;


if( nargin>=1 )
  arg1 = varargin{1};
  isVerbose = true;
  if( ischar(arg1) && length(arg1)>1 )
    sTest = arg1;
  end

  if( nargin>=2 )
    arg2 = varargin{2};
    if( ischar(arg2) )
      fid = fopen( arg2, 'a+' );
    end
  end
end
addpath( genpath('test/mloct_test') );


results = mloct_test_run( sSrcDir, sTestDir, sTest, fid, isVerbose );


if( fid>1 )
  fclose( fid );
end
if( ~nargout )
  clear results
end
