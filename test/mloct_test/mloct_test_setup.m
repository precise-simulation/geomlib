function [ fid, isVerbose, nCnt ] = mloct_test_setup( fid, isVerbose, nCnt )
%MLOCT_TEST_SETUP Test setup
%
% [ FID, ISVERBOSE, NCNT ] = MLOCT_TEST_SETUP( FID, ISVERBOSE, NCNT )
% Store, cache, and retrieve FID, ISVERBOSE, and NCNT test parameters.

% Initial version 180216.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.

persistent cFid cVerbose cCnt


if( isempty(cFid) )
  cFid = 1;
end
if( isempty(cVerbose) )
  cVerbose = false;
end
if( isempty(cCnt) )
  cCnt = 0;
end


if( nargin>=1 && isscalar(fid) && isnumeric(fid) )
  cFid = fid;
end
if( nargin>=2 && isscalar(isVerbose) )
  cVerbose = isVerbose;
end
if( nargin>=3 && isscalar(nCnt) )
  cCnt = nCnt;
end


if( nargout>=1 )
  fid = cFid;
end
if( nargout>=2 )
  isVerbose = cVerbose;
end
if( nargout>=3 )
  nCnt = cCnt;
end
