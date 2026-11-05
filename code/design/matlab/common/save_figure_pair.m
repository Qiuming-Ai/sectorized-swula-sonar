function save_figure_pair(fig, output_stem)
%SAVE_FIGURE_PAIR Save editable FIG and publication-ready 300-dpi PNG.
folder = fileparts(output_stem);
if ~isempty(folder) && ~exist(folder,'dir'), mkdir(folder); end
savefig(fig, output_stem + ".fig");
exportgraphics(fig, output_stem + ".png", 'Resolution', 300);
end
