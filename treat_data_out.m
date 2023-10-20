function [data_out, data_best] = treat_data_out (data_in)

data_out = struct2table(data_in); 
data_out.loss = round(data_out.loss, 3);
idx = ismember(data_out.loss, min(data_out.loss));
data_best = data_out(idx,:);
data_out = sortrows(data_out,'loss');