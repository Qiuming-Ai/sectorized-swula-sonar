function masks = sector_masks(cfg, angles_deg)
%SECTOR_MASKS Logical angle masks for the six stitched sectors.
n_sector = size(cfg.sector_edges_deg, 1);
masks = false(n_sector, numel(angles_deg));
for k = 1:n_sector
    lo = cfg.sector_edges_deg(k,1);
    hi = cfg.sector_edges_deg(k,2);
    if k < n_sector
        masks(k,:) = angles_deg >= lo & angles_deg < hi;
    else
        masks(k,:) = angles_deg >= lo & angles_deg <= hi;
    end
end
end
