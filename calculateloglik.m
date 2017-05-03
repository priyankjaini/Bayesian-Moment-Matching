function loglikelihoodData = calculateloglik(theta, phi, lambda, data)

fZero = theta*phi + (1-theta)*lambda;
fOne = 1 - fZero;

ind0 = find(data);

loglikelihoodData = (log(fZero)*length(ind0) + log(fOne)*(length(data) - length(ind0)))/length(data);
% loglikelihoodData = log(likelihood);