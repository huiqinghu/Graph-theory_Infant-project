function [mean_pos, mean_neg, mean_abs] = calculate_mean_FC(mat_data)

mat_data_pos = mat_data(~isnan(mat_data) & mat_data > 0);
mat_data_neg = mat_data(~isnan(mat_data) & mat_data < 0);
mat_data_abs = abs(mat_data); mat_data_abs = mat_data_abs(~isnan(mat_data_abs));
mean_pos = mean(mat_data_pos);
mean_neg = mean(mat_data_neg);
mean_abs = mean(mat_data_abs);
end

%%%%Written by Peter Coppola
