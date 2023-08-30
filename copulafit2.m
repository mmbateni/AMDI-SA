function [alphaHat,rmse] = copulafit2(u)
%function [alphaHat,rmse,alphaHat_bnd] = copulafit2(u)
%COPULAFIT Fit a parametric copula to data.
%
%   PARAMHAT = COPULAFIT(FAMILY, U) returns an estimate PARAMHAT of the copula
%   parameter for an Archimedean copula specified by FAMILY, given data in U.  U
%   is an N-by-2 matrix of values in (0,1), representing N points in the unit
%   square.  FAMILY is 'Clayton', 'Frank', or 'Gumbel'.
%
%   [...] = COPULAFIT(..., 'Alpha', ALPHA) returns an approximate 100(1-ALPHA)%
%   confidence interval for the parameter estimate.
%
%   COPULAFIT uses 'minimum squared error' to fit the copula to U.  When U contains
%   data transformed to the unit hypercube by parametric or nonparaetric estimates of their
%   marginal cumulative distribution functions. 

%
%   See also COPULAFIT, ECDF, COPULACDF, COPULAPDF, COPULARND.

%   References:
%      [1] Bouye, E., Durrleman, V., Nikeghbali, A., Riboulet, G., Roncalli, T.
%          (2000), "Copulas for Finance: A Reading Guide and Some Applications",
%          Working Paper, Groupe de Recherche Operationnelle, Credit Lyonnais.


[~,d] = size(u);

if ndims(u)~=2 || d<2
    error(message('stats:copulafit:InvalidDataDimensions'));
elseif ~all(all(0 < u & u < 1))
    error(message('stats:copulafit:DataOutOfRange'));
end

families = {'clayton','frank','gumbel'};
options=[];
methods = { 'MSE'};
method = { 'MSE'};
method = internal.stats.getParamVal(method,methods,'METHOD');
options = statset(statset('copulafit'), options);
    %% Calculating Empirical Copula for GOFs; Empirical calculated based on
    %  Aghakouchk's work on Non-Parametric Approach
    cop_empri=empcop(u);
    %----------------------------------------------------------------------
    %% Calculating Empirical Kendall Function Based on Nelson's Work
    K_c=empkend(u,cop_empri);
    for fa=1:3
 FAMILY=   families{fa};
 family = internal.stats.getParamVal(FAMILY,families,'FAMILY');
% case {'clayton' 'frank' 'gumbel'}
%     if nargout > 2
%         error(message('stats:copulafit:TooManyOutputs'));
%     elseif d > 2
%         error(message('stats:copulafit:TooManyDimensions'));
%     end

    switch family
    case 'clayton' % a.k.a. Cook-Johnson
        nloglf = @(alpha)rmse_clayton(alpha,cop_empri,K_c);
        lowerBnd = options.TolBnd;
    case 'frank'
        nloglf =@(alpha) rmse_frank(alpha,cop_empri,K_c);
        [~,lowerBnd] = bracket1D(nloglf,5,-5); % 'lower', search descending from -5
        if ~isfinite(lowerBnd)
            error(message('stats:copulafit:NoLowerBnd'));
        end
    case 'gumbel' % a.k.a. Gumbel-Hougaard
        nloglf = @(alpha) rmse_gumbel(alpha,cop_empri,K_c);
        lowerBnd = 1 + options.TolBnd;
    end
    [lowerBnd,upperBnd] = bracket1D(nloglf,lowerBnd,5); % 'upper', search ascending from 5
    if ~isfinite(upperBnd)
        error(message('stats:copulafit:NoUpperBnd'));
    end
    opts = optimset(options);
    [alphaHat(fa),~,err,output] = fminbnd(nloglf, lowerBnd, upperBnd, opts);
    %[alphaHat_bnd(fa),~,~,~] = fminbnd(nloglf, -50, 50, opts);
    if (err == 0)
        % fminbnd may print its own output text; in any case give something
        % more statistical here, controllable via warning IDs.
        if output.funcCount >= options.MaxFunEvals
            warning(message('stats:copulafit:EvalLimit'));
        else
            warning(message('stats:copulafit:IterLimit'));
        end
    elseif (err < 0)
        error(message('stats:copulafit:NoSolution'));
    end
    varargout{1} = alphaHat;
end
    rmse=[rmse_clayton(alphaHat(1),cop_empri,K_c),rmse_frank(alphaHat(2),cop_empri,K_c),rmse_gumbel(alphaHat(3),cop_empri,K_c)];
    varargout{2}=rmse;
    %varargout{3}=alphaHat_bnd;

%
% RMSE for Archimedean copulas 
%
    function [nll] = rmse_clayton(alpha,cop_empri,K_c)
    % C(u1,u2) = (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)
%           y_clayton =copulacdf('Clayton',u,alpha);
            y_clayton =cop_empri;
            K_c_clayton=y_clayton.*((1+alpha-(y_clayton.^alpha))/alpha);
            nll2=(sqrt((K_c(:)-K_c_clayton(:)).^2));
            nll1 = sum(nll2)./sqrt(numel(K_c));
            nll=nll1;
    end

    function [nll] = rmse_frank(alpha,cop_empri,K_c)
    % C(u1,u2) = -(1/alpha)*log(1 + (exp(-alpha*u1)-1)*(exp(-alpha*u1)-1)/(exp(-alpha)-1))
%               y_frank =  copulacdf('Frank',u,alpha);
                 y_frank =cop_empri;
                ex_K_frank=exp(-alpha*y_frank);
                K_c_frank=y_frank+(((1-ex_K_frank)./(alpha*ex_K_frank)).*(log((1-exp(-alpha))./(1-ex_K_frank))));
                nll2=sqrt((K_c(:)-K_c_frank(:)).^2) ;
                nll1 = sum(nll2)./sqrt(numel(K_c));
               nll=nll1;
    end

    function [nll] = rmse_gumbel(alpha,cop_empri,K_c)
    % C(u1,u2) = exp(-((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha))
%       y_gumbel = copulacdf('Gumbel',u,alpha);
        y_gumbel =cop_empri;
        K_c_gumbel=y_gumbel-((y_gumbel.*log(y_gumbel))/(alpha+1));
        nll2=sqrt((K_c(:)-K_c_gumbel(:)).^2);
        nll1 = sum(nll2)./sqrt(numel(K_c));
        nll=nll1;

    end

%----------------------------------------------------------------------




%----------------------------------------------------------------------
%
% Utility functions
%

% function x = tinvLocal(p,nu)
% 
% % For small d.f., call betaincinv which uses Newton's method
% if nu < 1000
%     q = p - .5;
%     z = zeros(size(q),class(q));
%     oneminusz = zeros(size(q),class(q));
%     t = (abs(q(:)) < .25);
%     if any(t)
%         % for z close to 1, compute 1-z directly to avoid roundoff
%         oneminusz(t) = betaincinv(2.*abs(q(t)),0.5,nu/2,'lower');
%         z(t) = 1 - oneminusz(t);
%     end
%     t = ~t; % (abs(q) >= .25);
%     if any(t)
%         z(t) = betaincinv(2.*abs(q(t)),nu/2,0.5,'upper');
%         oneminusz(t) = 1 - z(t);
%     end
%     x = sign(q) .* sqrt(nu .* (oneminusz./z));
%     
% % For large d.f., use Abramowitz & Stegun formula 26.7.5
% else
%     xn = norminv(p);
%     x = xn + (xn.^3+xn)./(4*nu) + ...
%         (5*xn.^5+16.*xn.^3+3*xn)./(96*nu.^2) + ...
%         (3*xn.^7+19*xn.^5+17*xn.^3-15*xn)./(384*nu.^3) +...
%         (79*xn.^9+776*xn.^7+1482*xn.^5-1920*xn.^3-945*xn)./(92160*nu.^4);
% end
% end

function [nearBnd,farBnd] = bracket1D(nllFun,nearBnd,farStart)
% Bracket the minimizer of a (one-param) negative log-likelihood function.
% nearBnd is a point known to be a lower/upper bound for the minimizer,
% this will be updated to tighten the bound if possible.  farStart is the
% first trial point to test to see if it's an upper/lower bound for the
% minimizer.  farBnd will be the desired upper/lower bound.
bound = farStart;
upperLim = 1e12; % arbitrary finite limit for search
oldnll = nllFun(bound);
oldbound = bound;
while abs(bound) <= upperLim
    bound = 2*bound; % assumes lower start is < 0, upper is > 0
    nll = nllFun(bound);
    if nll > oldnll
        % The neg loglikelihood increased, we're on the far side of the
        % minimum, so the current point is the desired far bound.
        farBnd = bound;
        break;
    else
        % The neg loglikelihood continued to decrease, so the previous point
        % is on the near side of the minimum, update the near bound.
        nearBnd = oldbound;
    end
    oldnll = nll;
    oldbound = bound;
end
if abs(bound) > upperLim
    farBnd = NaN;
end
end


function b = tovector(A)
% Convert the Cholesky factor of a correlation matrix to upper triangle vector
% form.  Columns of the Cholesky factor have unit norm, so they reduce to
% direction angles, consisting of one fewer element.
m = size(A,1);
Angles = zeros(m);
for i = 1:m-1
    Angles(i,:) = atan2(A(i,:),A(i+1,:));
    A(i+1,:) = A(i+1,:) ./ cos(Angles(i,:));
end
b = Angles(triu(true(m),1));
end


function A = tomatrix(b)
% Convert the Cholesky factor of a correlation matrix from upper triangle vector
% form.  Columns of the Cholesky factor have unit norm, so the columns can be
% recreated from direction angles.
m = (1 + sqrt(1+8*length(b)))/2;
Cosines = zeros(m);
Cosines(triu(true(m),1)) = cos(b);
Sines = ones(m);
Sines(triu(true(m),1)) = sin(b);
prodSines = cumprod(Sines, 1, 'reverse'); 
A = [ones(1,m); Cosines(1:m-1,:)] .* prodSines;
% A = [ones(1,m); Cosines(1:m-1,:)] .* flipud(cumprod(flipud(Sines)));
end
end