function paramNew = CopyDefaultStructParams(paramDefault, param)
% Copy default parameter struct into a given user parameter struct
%
% Copyright (c) 2014, Jian Cheng (jian.cheng.1983@gmail.com)
%

%%
fields = fieldnames(paramDefault);
paramNew = param;

for i = 1 : numel(fields)
    if ~isfield(param, fields{i})
        paramNew = setfield(paramNew, fields{i}, getfield(paramDefault, fields{i}));
    end
end

