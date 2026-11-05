function [data, count] = legacy_load_timing_fixture(file, expected_type, max_pings)
% Input loading and checking occur before the timing loop.
assert(isfile(file), ['Missing timing HDF5 input: %s. Extract the optional ' ...
    'desktop-timing data overlay or pass the exact input file explicitly.'], file);
assert(isnumeric(max_pings) && isscalar(max_pings) && max_pings > 0 && ...
    (isinf(max_pings) || max_pings == floor(max_pings)), 'max_pings must be a positive integer or Inf.');
indices = [];
if ~isinf(max_pings)
    indices = 1:max_pings;
end
data = ReadBaselineHDF5(file, 'PingIndices', indices);
assert(isfield(data.attributes, 'array_type') && ...
    strcmpi(string(data.attributes.array_type), expected_type), 'Unexpected array type in %s.', file);
assert(isfield(data.attributes, 'ping_num') && isfield(data.attributes, 'receive_array_num'), ...
    'Incomplete timing input attributes in %s.', file);
count = min(double(data.attributes.ping_num), max_pings);
assert(isscalar(count) && isfinite(count) && count > 0 && count == floor(count) && ...
    numel(data.pings) == count, 'Missing or incomplete ping data in %s.', file);
assert(all(cellfun(@(x) size(x,2) == data.attributes.receive_array_num, data.pings)), ...
    'Ping channel dimensions disagree with input attributes: %s.', file);
end
