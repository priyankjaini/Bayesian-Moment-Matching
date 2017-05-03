function [alpha, beta, gamma] = momentMatchingSufficient(alpha, beta, gamma, data)

momentTheta(1) = alpha(1)/sum(alpha); 
momentTheta(2) = alpha(1).*(alpha(1)+1)/(sum(alpha)*(1+sum(alpha)));
momentTheta(3) = alpha(1).*(alpha(1)+1).*(alpha(1)+2)/(sum(alpha)*(1+sum(alpha))*(2+sum(alpha)));

momentPhi(1) = beta(1)/sum(beta); 
momentPhi(2) = beta(1).*(beta(1)+1)/(sum(beta)*(1+sum(beta)));
momentPhi(3) = beta(1).*(beta(1)+1).*(beta(1)+2)/(sum(beta)*(1+sum(beta))*(2+sum(beta)));

momentLambda(1) = gamma(1)/sum(gamma); 
momentLambda(2) = gamma(1).*(gamma(1)+1)/(sum(gamma)*(1+sum(gamma)));
momentLambda(3) = gamma(1).*(gamma(1)+1).*(gamma(1)+2)/(sum(gamma)*(1+sum(gamma))*(2+sum(gamma)));

if (data == 0)
    %%% Normalization Constant P(F=0) %%%
    Z = momentTheta(1)*momentPhi(1) + momentLambda(1) - momentTheta(1)*momentLambda(1);
    %%% Set of Minimal Moments %%%
    lambda(1,1) = (momentTheta(1)*momentPhi(1)*momentLambda(1) + momentLambda(2) - momentTheta(1)*momentLambda(2))/Z;
    theta(1,1) = (momentTheta(2)*momentPhi(1) + momentTheta(1)*momentLambda(1) - momentTheta(2)*momentLambda(1) )/Z;
    phi(1,1) = ( momentTheta(1)*momentPhi(2) + momentPhi(1)*momentLambda(1) - momentTheta(1)*momentPhi(1)*momentLambda(1) )/Z;
    %%% Check if the moments preserve the property E(XY) < E(X) and < E(Y) %%%
    
    %%% Extra Moments %%%
    theta(1,2) = (momentTheta(3)*momentPhi(1) + momentTheta(2)*momentLambda(1) - momentTheta(3)*momentLambda(1))/Z;
    phi(1,2) = (momentTheta(1)*momentPhi(3) + momentPhi(2)*momentLambda(1) - momentTheta(1)*momentPhi(2)*momentLambda(1))/Z;
    lambda(1,2) = (momentTheta(1)*momentPhi(1)*momentLambda(2) + momentLambda(3) - momentTheta(1)*momentLambda(3))/Z;
    
    %%%
    theta(2,1) = 1 - theta(1,1);
    phi(2,1) = 1 - phi(1,1);
    lambda(2,1) = 1 - lambda(1,1);
    theta(2,2) = 1 + theta(1,2) - 2*theta(1,1);
    phi(2,2) = 1 + phi(1,2) - 2*phi(1,1);
    lambda(2,2) = 1 + lambda(1,2) - 2*lambda(1,1);
    %%% Projection %%%
    alpha = theta(:,1).*(theta(:,1) - theta(:,2))./(theta(:,2) - theta(:,1).^2);
    beta = phi(:,1).*(phi(:,1) - phi(:,2))./(phi(:,2) - phi(:,1).^2);
    gamma = lambda(:,1).*(lambda(:,1) - lambda(:,2))./(lambda(:,2) - lambda(:,1).^2);
else
    %%% Normalization Constant P(F=1) %%%
    Z = 1 - (momentTheta(1)*momentPhi(1) + momentLambda(1) - momentTheta(1)*momentLambda(1));
    %%% Set of Sufficient Moments %%%
    lambda(1,1) = ( momentLambda(1) - (momentTheta(1)*momentPhi(1)*momentLambda(1) + momentLambda(2) - momentTheta(1)*momentLambda(2)) )/Z;
    theta(1,1) = ( momentTheta(1) - (momentTheta(2)*momentPhi(1) + momentTheta(1)*momentLambda(1) - momentTheta(2)*momentLambda(1)) )/Z;
    phi(1,1) = ( momentPhi(1) - (momentTheta(1)*momentPhi(2) + momentPhi(1)*momentLambda(1) - momentTheta(1)*momentPhi(1)*momentLambda(1)) )/Z;
    %%% Check if the moments preserve the property E(XY) < E(X) and < E(Y) %%%
    
    %%% Extra Moments %%%
    theta(1,2) = (momentTheta(2) - (momentTheta(3)*momentPhi(1) + momentTheta(2)*momentLambda(1) - momentTheta(3)*momentLambda(1)) )/Z;
    phi(1,2) = (momentPhi(2) - (momentTheta(1)*momentPhi(3) + momentPhi(2)*momentLambda(1) - momentTheta(1)*momentPhi(2)*momentLambda(1)) )/Z;
    lambda(1,2) = (momentLambda(2) - (momentTheta(1)*momentPhi(1)*momentLambda(2) + momentLambda(3) - momentTheta(1)*momentLambda(3)) )/Z;
    
    %%%
    theta(2,1) = 1 - theta(1,1);
    phi(2,1) = 1 - phi(1,1);
    lambda(2,1) = 1 - lambda(1,1);
    theta(2,2) = 1 + theta(1,2) - 2*theta(1,1);
    phi(2,2) = 1 + phi(1,2) - 2*phi(1,1);
    lambda(2,2) = 1 + lambda(1,2) - 2*lambda(1,1);
    %%% Projection %%%
    alpha = theta(:,1).*(theta(:,1) - theta(:,2))./(theta(:,2) - theta(:,1).^2);
    beta = phi(:,1).*(phi(:,1) - phi(:,2))./(phi(:,2) - phi(:,1).^2);
    gamma = lambda(:,1).*(lambda(:,1) - lambda(:,2))./(lambda(:,2) - lambda(:,1).^2);
    
end