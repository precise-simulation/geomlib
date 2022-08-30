function [ subFcns ] = get_subfuns( sFun, sPrefix )
%GET_SUBFUNS Get sub function names
%
% [ SUBFCNS ] = GET_SUBFUNS( SFUN, SPREFIX ) Retreive subfunction
% names in SFUN (optionally filter out subfunctions without SPREFIX).

% Initial version 180216.
% Copyright 2013-2022 Precise Simulation Ltd.
% License: AGPL v3, see LICENSE for more details or contact
%          Precise Simulation for alternative licensing options.


sFile = fileread( sFun );
cLines = regexp( strtrim(sFile), '\n', 'split' );

cLines( ~strncmpi(cLines,'function',8) ) = [];
cLines(1) = [];


for i=1:length(cLines)
  sSubFun = strtrim(cLines{i});
  i1 = max([9,find(sSubFun== '=',1)+1]);
  i2 = min([length(sSubFun),find(sSubFun== '(',1)-1]);
  cLines{i} = strtrim(sSubFun(i1:i2));
end


subFcns = cLines;
if( nargin>=2 )
  for i=length(subFcns):-1:1
    if( ~strncmpi(subFcns{i},'test',4) )
      subFcns(i) = [];
    end
  end
end
