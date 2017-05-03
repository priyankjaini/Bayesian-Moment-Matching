function [alpha, beta, gamma] = momentMatchingMartingale(alpha, beta, gamma, data)


%%% FIRST ORDER MOMENTS %%%
momentTheta(1) = alpha(1)/sum(alpha); 
momentTheta(2) = alpha(1).*(alpha(1)+1)/(sum(alpha)*(1+sum(alpha)));
momentTheta(3) = alpha(1).*(alpha(1)+1).*(alpha(1)+2)/(sum(alpha)*(1+sum(alpha))*(2+sum(alpha)));

momentPhi(1) = beta(1)/sum(beta); 
momentPhi(2) = beta(1).*(beta(1)+1)/(sum(beta)*(1+sum(beta)));
momentPhi(3) = beta(1).*(beta(1)+1).*(beta(1)+2)/(sum(beta)*(1+sum(beta))*(2+sum(beta)));

momentLambda(1) = gamma(1)/sum(gamma); 
momentLambda(2) = gamma(1).*(gamma(1)+1)/(sum(gamma)*(1+sum(gamma)));
momentLambda(3) = gamma(1).*(gamma(1)+1).*(gamma(1)+2)/(sum(gamma)*(1+sum(gamma))*(2+sum(gamma)));

% T = alpha/sum(alpha);
% P = beta/sum(beta);
% L = gamma/sum(gamma);
% 
% %%% SECOND ORDER MOMENTS %%%
% T2 = alpha.*(alpha+1)/(sum(alpha)*(1+sum(alpha)));
% P2 = beta.*(beta+1)/(sum(beta)*(1+sum(beta)));
% L2 = gamma.*(gamma+1)/(sum(gamma)*(1+sum(gamma)));
% 
% %%% THIRD ORDER MOMENTS %%%
% T3 = alpha.*(alpha+1).*(alpha+2)/(sum(alpha)*(1+sum(alpha))*(2+sum(alpha)));
% P3 = beta.*(beta+1).*(beta+2)/(sum(beta)*(1+sum(beta))*(2+sum(beta)));
% L3 = gamma.*(gamma+1).*(gamma+2)/(sum(gamma)*(1+sum(gamma))*(2+sum(gamma)));

if (data == 0)
    %%% Normalization Constant P(F=0) %%%
    Z = momentTheta(1)*momentPhi(1) + momentLambda(1) - momentTheta(1)*momentLambda(1);
    %%% Set of Minimal Moments %%%
    theta_phi = (momentTheta(2)*momentPhi(2) + momentTheta(1)*momentPhi(1)*momentLambda(1) - momentTheta(2)*momentPhi(1)*momentLambda(1))/Z;
    theta_lambda = (momentTheta(2)*momentPhi(1)*momentLambda(1) + momentTheta(1)*momentLambda(2) - momentTheta(2)*momentLambda(2))/Z;
    lambda(1,1) = (momentTheta(1)*momentPhi(1)*momentLambda(1) + momentLambda(2) - momentTheta(1)*momentLambda(2))/Z;
    
    theta(1,1) = theta_lambda/lambda(1);
    phi(1,1) = theta_phi/theta(1);
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
    %%% Set of Minimal Moments %%%
    theta_phi = ( momentTheta(1)*momentPhi(1) - (momentTheta(2)*momentPhi(2) + momentTheta(1)*momentPhi(1)*momentLambda(1) - momentTheta(2)*momentPhi(1)*momentLambda(1)) )/Z;
    theta_lambda = ( momentTheta(1)*momentLambda(1) - (momentTheta(2)*momentPhi(1)*momentLambda(1) + momentTheta(1)*momentLambda(2) - momentTheta(2)*momentLambda(2)) )/Z;
    lambda(1,1) = ( momentLambda(1) - (momentTheta(1)*momentPhi(1)*momentLambda(1) + momentLambda(2) - momentTheta(1)*momentLambda(2)) )/Z;
    
    theta(1,1) = theta_lambda/lambda(1);
    phi(1,1) = theta_phi/theta(1);
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
    
    %         theta(1,1) = (T2(1)*P(1) + alpha(1)*alpha(2)*L(1)/( sum(alpha)*(sum(alpha)+1) ))/Z;
    %         theta(2,1) = (T2(2)*L(1) + alpha(1)*alpha(2)*P(1)/( sum(alpha)*(sum(alpha)+1) ))/Z;
    %         %%% MOMENTS TO MATCH FOR A MARTINGALE : { theta1_phi1, theta2_lambda1 }
    %
    %         theta1_phi1 = ( T2(1)*momentPhi(2) + ( alpha(1)*alpha(2)*P(1)*L(1) )/( sum(alpha)*(1+sum(alpha)) ) )/Z;
    %         theta2_lambda1 = ( T2(2)*momentLambda(2) + ( alpha(1)*alpha(2)*P(1)*L(1) )/( sum(alpha)*(1+sum(alpha)) ) )/Z;
    %
    %         %%% REST OF THE SUFFICIENT MOMENTS
    %                 phi(1,1) = T(1)*momentPhi(2) + T(2)*P(1)*L(1);
    % %%%         phi(1,1) = theta1_phi1/theta(1,1);
    %         phi(1,2) = T(1)*momentPhi(3) + T(2)*momentPhi(2)*L(1);
    %         phi(2,1) = ( T(1)*beta(1)*beta(2) )/( sum(beta)*(1+sum(beta)) ) + T(2)*L(1)*P(2);
    %         phi(2,2) = ( T(1)*beta(1)*beta(2)*(beta(2)+1) )/( sum(beta)*(1+sum(beta))*(2+sum(beta)) ) + T(2)*L(1)*(P2(2));
    %         phi = phi/Z;
    %
    %                 lambda(1,1) = T(1)*P(1)*L(1) + T(2)*momentLambda(2);
    % %%%         lambda(1,1) = theta2_lambda1/theta(2,1);
    %         lambda(1,2) = T(1)*P(1)*momentLambda(2) + T(2)*momentLambda(3);
    %         lambda(2,1) = T(1)*P(1)*L(2) + ( T(2)*gamma(1)*gamma(2) )/( sum(gamma)*(1+sum(gamma)) );
    %         lambda(2,2) = T(1)*P(1)*L2(2) + ( T(2)*gamma(1)*gamma(2)*(gamma(2)+1) )/( sum(gamma)*(1+sum(gamma))*(2+sum(gamma)) );
    %         lambda = lambda/Z;
    %
    %         %%%         theta(1,1) = theta1_phi1/phi(1,1);
    %         %%%         theta(2,1) = theta2_lambda1/lambda(1,1);
    %
    %         theta(1,2) = ( momentTheta(3)*P(1) + L(1)*alpha(1)*alpha(2)*(alpha(1)+1)/(sum(alpha)*(1+sum(alpha))*(2+sum(alpha))) )/Z;
    %         theta(2,2) = ( ( P(1)*alpha(1)*alpha(2)*(alpha(2)+1) )/( sum(alpha)*(1+sum(alpha))*(2+sum(alpha)) ) + T3(2)*L(1) )/Z;
    %         %%%         theta = theta/Z;
    %         %%%% PROJECTION %%%%
    %         alpha = theta(:,1).*(theta(:,1) - theta(:,2))./( theta(:,2) - theta(:,1).^2 );
    %         beta = phi(:,1).*(phi(:,1) - phi(:,2))./( phi(:,2) - phi(:,1).^2 );
    %         gamma = lambda(:,1).*(lambda(:,1) - lambda(:,2))./( lambda(:,2) - lambda(:,1).^2 );
    %     else
    %         Z = T(1)*P(2) + T(2)*L(2);
    %         theta(1,1) = ( T2(1)*P(2) + alpha(1)*alpha(2)*L(2)/( sum(alpha)*(sum(alpha)+1) ))/Z;
    %         theta(2,1) = (T2(2)*L(2) + alpha(1)*alpha(2)*P(2)/( sum(alpha)*(sum(alpha)+1) ))/Z;
    %         %%% MOMENTS TO MATCH FOR A MARTINGALE : { theta1_phi2, theta2_lambda2 }
    %         theta1_phi2 = ( T2(1)*P2(2) + ( alpha(1)*alpha(2)*P(2)*L(2) )/( sum(alpha)*(1+sum(alpha)) ) )/Z;
    %         theta2_lambda2 = ( T2(2)*L2(2) + ( alpha(1)*alpha(2)*P(2)*L(2) )/( sum(alpha)*(1+sum(alpha)) ) )/Z;
    %         %%% REST OF THE SUFFICIENT MOMENTS
    %         phi(1,1) = ( T(1)*beta(1)*beta(2) )/( sum(beta)*(1+sum(beta)) ) + T(2)*P(1)*L(2);
    %         phi(1,2) = ( T(1)*beta(1)*beta(2)*(beta(1)+1) )/( sum(beta)*(1+sum(beta))*(2+sum(beta)) ) + T(2)*momentPhi(2)*L(2);
    %         phi(2,1) = T(1)*P2(2) + T(2)*L(2)*P(2);
    % %%%%         phi(2,1) = theta1_phi2/theta(1,1);
    %         phi(2,2) = T(1)*P3(2) + T(2)*L(2)*P2(2);
    %         phi = phi/Z;
    %
    %         lambda(1,1) = T(1)*P(2)*L(1) + ( T(2)*gamma(1)*gamma(2) )/( sum(gamma)*(1+sum(gamma)) );
    %         lambda(1,2) = T(1)*P(2)*momentLambda(2) + ( T(2)*gamma(1)*gamma(2)*(gamma(1)+1) )/( sum(gamma)*(1+sum(gamma))*(2+sum(gamma)) );
    %         lambda(2,1) = T(1)*P(2)*L(2) + T(2)*L2(2);
    % %%%         lambda(2,1) = theta2_lambda2/theta(2,1);
    %         lambda(2,2) = T(1)*P(2)*L2(2) + T(2)*L3(2);
    %         lambda = lambda/Z;
    %
    %         %%%         theta(1,1) = theta1_phi2/phi(2,1);
    %         %%%         theta(2,1) = theta2_lambda2/lambda(2,1);
    %
    %         theta(1,2) = ( momentTheta(3)*P(2) + ( L(2)*alpha(1)*alpha(2)*(alpha(1)+1) )/( sum(alpha)*(1+sum(alpha))*(2+sum(alpha)) ) )/Z;
    %         theta(2,2) = ( ( P(2)*alpha(1)*alpha(2)*(alpha(2)+1) )/( sum(alpha)*(1+sum(alpha))*(2+sum(alpha)) ) + T3(2)*L(2) )/Z;
    %%%% PROJECTION %%%%
    %         alpha = theta(:,1).*(theta(:,1) - theta(:,2))./( theta(:,2) - theta(:,1).^2 );
    %         beta = phi(:,1).*(phi(:,1) - phi(:,2))./( phi(:,2) - phi(:,1).^2 );
    %         gamma = lambda(:,1).*(lambda(:,1) - lambda(:,2))./( lambda(:,2) - lambda(:,1).^2 );
end