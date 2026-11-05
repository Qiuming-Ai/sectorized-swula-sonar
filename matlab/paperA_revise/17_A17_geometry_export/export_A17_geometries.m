function geometry_table=export_A17_geometries()
%EXPORT_A17_GEOMETRIES Exact legacy sparse-array coordinates used in Fig. 16.
this_dir=fileparts(mfilename('fullpath'));
out_dir=fullfile(this_dir,'results');if ~exist(out_dir,'dir'),mkdir(out_dir);end
c0=1500;f0=450e3;lambda=c0/f0;pitch=0.5*lambda;

% The figure export injects rng(1401,'twister') before the legacy script.
rng(1401,'twister');n0=128;n_random=44;
candidate=2:n0-1; pick=candidate(randperm(numel(candidate),n_random-2));
all_dense=linspace(-n0*pitch/2,n0*pitch/2,n0).';
random_m=all_dense([1,pick,128]);

n1=12;n2=32;pna=n1*n2;d0=63.5*lambda/pna;
nested_m=[(0:n1-1)'*d0;(1:n2)'*n1*d0-63.5*lambda/2];

n1c=17;n2c=27;pnac=n1c*n2c;d0c=63.5*lambda/pnac;
coprime_m=[(0:n1c-1)'*n2c*d0c-63.5*lambda/2;(1:n2c)'*n1c*d0c-63.5*lambda/2];

method=[repmat("random",numel(random_m),1);repmat("nested",numel(nested_m),1);repmat("coprime",numel(coprime_m),1)];
element_index=[(1:numel(random_m))';(1:numel(nested_m))';(1:numel(coprime_m))'];
position_m=[random_m;nested_m;coprime_m];position_lambda=position_m/lambda;
geometry_table=table(method,element_index,position_m,position_lambda);
writetable(geometry_table,fullfile(out_dir,'A17_exact_geometry_coordinates.csv'));
save(fullfile(out_dir,'A17_geometry_workspace.mat'),'geometry_table');
disp(groupsummary(geometry_table,'method',{'min','max'},'position_m'));
end
