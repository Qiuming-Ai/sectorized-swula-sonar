function cfg=config_A16(mode,geometry_dir)
%CONFIG_A16 Phase-ensemble audit using the saved A11 arrival cloud.
if nargin<1,mode="quick";end
if nargin<2,geometry_dir=[];end
mode=validatestring(mode,{'quick','paper'});
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'common'));
a11_dir=fullfile(root,'11_A11_shipwreck_common_scene');addpath(a11_dir);
cfg=config_A11("paper",geometry_dir,"archive");
cfg.methods=["dense","cdm","fdm_current"];
cfg.method_labels=["Dense ULA","CDM SW-ULA","packed FDM SW-ULA"];
if lower(string(mode))=="quick",cfg.phase_runs=3;else,cfg.phase_runs=30;end
cfg.phase_seed_base=cfg.ship_seed;
cfg.noise_seed_base=2026081600;
cfg.guard_radii=[4,1;8,1;16,2;32,4];
cfg.primary_guard_radius=[16,2];
package_root=fileparts(fileparts(root));
cfg.output_dir=fullfile(package_root,'outputs','matlab',"A16_"+string(mode));
cfg.saved_a11_workspace="";
cfg.phase_definition="independent uniform unresolved-scatterer phase on each saved one-bounce arrival; phases are shared by all imaging methods within a run";
end
