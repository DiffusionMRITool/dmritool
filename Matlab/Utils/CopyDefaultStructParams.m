function paramReal = CopyDefaultStructParams(paramDefault, param)
% Copy default parameter struct into a given user parameter struct
%
% USAGE:
%    paramReal = CopyDefaultStructParams(paramDefault, param)
%
% INPUT
%    paraDefault  :  Parameter struct with all default values.
%    param        :  The input parameter  
%
% INPUT
%    paramReal    :  Real parameter struct used in the following codes with a combination of paramDefault and param.
%                    The default values in paramDefault are overriden by param
%
% Copyright (c) 2014, Jian Cheng (jian.cheng.1983@gmail.com)
%

%%
fields = fieldnames(paramDefault);
paramReal = param;

for i = 1 : numel(fields)
    if ~isfield(param, fields{i})
        paramReal = setfield(paramReal, fields{i}, getfield(paramDefault, fields{i}));
    end
end

