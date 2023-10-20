function [out_range_min, out_range_diap] = fit_range (in_range)

out_range_min  = min(in_range);
out_range_diap = max(in_range)-min(in_range);

