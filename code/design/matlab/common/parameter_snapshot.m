function parameter_snapshot(cfg, filename)
%PARAMETER_SNAPSHOT Save the exact configuration used by a run as JSON.
text = jsonencode(cfg, PrettyPrint=true);
fid = fopen(filename, 'w');
if fid < 0, error('Cannot write %s', filename); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);
end
