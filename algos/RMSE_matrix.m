function out = RMSE_matrix(X)

% X = X(~isnan(X));
for i =1:size(X,2)
    temp = X(~isnan(X(:,i)),i);
    out(i) = sqrt(mean(temp.^2));
end
