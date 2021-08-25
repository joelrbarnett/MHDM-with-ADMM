function w = ADMM_literature_tight(f, lam, alp, xk,w0)
% this function solves min_w lam* (sum_i wi + xk_i + f_i exp(-w_i-xk_i)) + lam*
% alp*TV(w+xk)+ TV(w)
% % INPUTS: 
%       f:      n1*n2 matrix, noisy image
%       lam:    positive scalar, multiscale parameter
%       alp:    positive scalar, the parameter in front of TV
%       w0:     n1*n2 matrix, initialization
% % OUTPUT:
%       w:      n1*n2 matrix, the log restored image piece
%
N = 100000; epsl = 1e-4; % stopping criteria
rho = 1.0;
y0 = 0.*f; theta0=w0; psi1_0=w0; psi2_0=w0;
for i=1:N
	% theta1 solves min_theta lam*(sum_i theta_i+ xk_i + f_i exp(-theta_i-xk_i)) + rho/2(theta - (psi1+pis2)/2+y)^2
	% by Newton's method solving rho*theta_i-lam*f_i*exp(-theta_i-xk_i) =
	% -lam + rho((psi1+psi2)_i/2-y_i)
	theta1 = Newton_step(theta0,f,xk,lam,-lam+rho*((psi1_0+psi2_0)./2-y0),rho, epsl * 0.1);
    
    %solve for minimizer, psi, of 
    % alp*lam/rho*TV(psi+xk) + 1/2*(theta1-(psi+psi2_0)/2+y0)^2
    psi1_1 = TVL2(2*theta1+2*y0+xk-psi2_0,4*alp*lam/rho , 4, 0)-xk;
    
    %find min of TV(psi) + rho/2*(theta1-(psi1_1+psi)/2+y0)^2 by finding
    %minimizer of 4TV(psi)/rho+ 1/2*(2*theta1-psi1_1+2*y0-psi)^2 via TVL2
    psi2_1 = TVL2(2*theta1-psi1_1+2*y0,4/rho,4,0);
    
    %update residual
    y1 = y0 + theta1 - (psi1_1+psi2_1)./2;
    
    %compute errors for exit condition
    err = max([norm(y1(:)-y0(:)); norm(theta1(:)-theta0(:)); norm(psi1_1(:)-psi1_0(:)); norm(psi2_1(:)-psi2_0(:)); norm(theta1(:)-(psi1_1(:)+psi2_1(:))./2)].^2);
    if err < epsl
        break;
    end
    fprintf('iter %d: error: %.4f\n', i, err);
    y0 = y1; theta0=theta1; psi1_0 = psi1_1; psi2_0=psi2_1;
end
w=theta1;
end

function p = Newton_step(p0, f, xk, lam, c, rho, epsl)
% this function implements Newton's method for rho*p-lam*f*exp(-p-xk) = c (pointwisely)
% % INPUTS: 
%       p0:     n1*n2 matrix, initialization
%       lam:    scalar
%       rho:    scalar
%       c:      n1*n2 matrix
%       epsl:   stopping tolerate
% % OUTPUT:
%       p:      n1*n2 matrix, solution
%
N = 100000;
for i=1:N
    p1 = p0 - (rho.*p0-lam*f.*exp(-p0-xk) - c) ./ (rho+lam*f.*exp(-p0-xk));
    err = sum((p1(:)-p0(:)).^2);
    if err < epsl
        %fprintf('error in Newtons is: %d . \n',err);
        break;
    end
    p0 = p1;
end
if i == N
    fprintf('error in Newtons is: %d . \n',err);
    fprintf('Newton does not converge.\n');
end
p = p1;
end