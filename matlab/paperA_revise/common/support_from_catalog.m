function support = support_from_catalog(catalog,n_range,n_angle)
%SUPPORT_FROM_CATALOG Logical truth occupancy from saved scatterer coordinates.
support = false(n_range,n_angle);
range_index = round(catalog.RangeIndex(:));
angle_index = round(catalog.AngleIndex(:));
valid = range_index>=1 & range_index<=n_range & ...
    angle_index>=1 & angle_index<=n_angle;
linear = sub2ind([n_range,n_angle],range_index(valid),angle_index(valid));
support(linear) = true;
end
