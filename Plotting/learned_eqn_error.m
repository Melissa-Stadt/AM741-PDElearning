function [err_array, err_mean] = learned_eqn_error(true_coeff, learned_coeff)
% Input: array of true_coeff and learned_coeff
% NOTE: if did not learn all coeff set 0 for those not learned
%   put coefficients that match the true_coeff first and then add
%   the other coefficients at the end
% Output: err array if matching
err_array = zeros(length(true_coeff), 1);
for ii = 1:length(true_coeff)
    err_temp=abs((true_coeff(ii)-learned_coeff(ii))/true_coeff(ii));
    err_array(ii) = err_temp;
end % for
if length(learned_coeff) > length(true_coeff)
    added_error = ones(length(learned_coeff) - length(true_coeff), 1);
    err_array = [err_array; added_error];
end % if

% mean of error
err_mean = mean(err_array);
end % end of learned eqn error