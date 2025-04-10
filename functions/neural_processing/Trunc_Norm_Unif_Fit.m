function [paramEsts, negloglik] = Trunc_Norm_Unif_Fit(x,K)

truncVar = 2;


% k means initialization.
if K == 1
    [cInd, C] = kmeans(x, K, 'EmptyAction','singleton');
    [Y,I] = sort(C);
else
    [cInd, C] = kmeans(x, K-1, 'EmptyAction','singleton');
    [Y,I] = sort(C);
end

% initial parameters for mixture fit. means and standard deviations
% from kmeans.
if K == 1
    start = [Y', std(x(cInd==I(1)))];
elseif K == 2
    start = [.5, Y', std(x(cInd==I(1)))]; % for 2
elseif K == 3
    start = [1/3, 1/3, Y', std(x(cInd==I(1))), std(x(cInd==I(2)))]; % if K  3
elseif K == 4
    start = [1/4, 1/4, 1/4, Y', std(x(cInd==I(1))), std(x(cInd==I(2))), std(x(cInd==I(3)))]; % if K  4
end


%% defining the type of distribution to fit

% normal distribution truncated on left and right side
xTrunc = [0 0.5];
switch truncVar
    case 0
        pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma);
    case 1
        pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma) ./ ...
            (1-normcdf(xTrunc(1),mu,sigma));
    case 2
        pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma) ./ ...
            (normcdf(xTrunc(2),mu,sigma)-normcdf(xTrunc(1),mu,sigma));
end

if K == 1
    % mixture of 1 truncated normals
    pdf_normmixture = @(x,mu1,sigma1) ...
        pdf_truncnorm(x,mu1,sigma1);
    
elseif K == 2
    % mixture of 2 truncated normals
    pdf_normmixture = @(x,p,mu1,sigma1) ...
        p*pdf_truncnorm(x,mu1,sigma1) + ...
        (1-p)*unifpdf(x,xTrunc(1),xTrunc(2));
    
elseif K == 3
    % mixture of 3 truncated normals
    pdf_normmixture = @(x,p1,p2,mu1,mu2,sigma1,sigma2) ...
        p1*pdf_truncnorm(x,mu1,sigma1) + ...
        p2*pdf_truncnorm(x,mu2,sigma2) + ...
        (1-(p1+p2))*unifpdf(x,xTrunc(1),xTrunc(2));
    
elseif K == 4
    % mixture of 3 truncated normals
    pdf_normmixture = @(x,p1,p2,p3,mu1,mu2,mu3,sigma1,sigma2,sigma3) ...
        p1*pdf_truncnorm(x,mu1,sigma1) + ...
        p2*pdf_truncnorm(x,mu2,sigma2) + ...
        p3*pdf_truncnorm(x,mu3,sigma3) + ...
        (1-(p1+p2+p3))*unifpdf(x,xTrunc(1),xTrunc(2));
end
% for troubleshooting
% t = 0:.005:.5;
% plot(t,pdf_normmixture(t,start(1),start(2),start(3),start(4),start(5),start(6),start(7),start(8)));


%% defining starting parameters for the maximum likelihood estimation

% pStart = .5;
% use kmeans here.
% start = [pStart mu1Start mu2Start sigma1Start sigma2Start];

lb = zeros(size(start));
if K == 1
    ub = [1 .5*ones(1,length(start)-1)];
else
    ub = [ones(1,K-1) .5*ones(1,length(start)-(K-1))];
end
options = statset('MaxIter',10000, 'MaxFunEvals',5000);
try
    paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, 'lower',lb,'upper',ub,'options',options);
catch
    paramEsts = nan;
    negloglik = nan;
    return
end

if K == 1
    negloglik = -sum(log(pdf_normmixture(x,paramEsts(1),paramEsts(2))));
elseif K == 2
    negloglik = -sum(log(pdf_normmixture(x,paramEsts(1),paramEsts(2),paramEsts(3))));
elseif K == 3
    negloglik = -sum(log(pdf_normmixture(x,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5),paramEsts(6))));
elseif K ==4
    negloglik = -sum(log(pdf_normmixture(x,paramEsts(1),paramEsts(2),paramEsts(3),...
        paramEsts(4),paramEsts(5),paramEsts(6),paramEsts(7),paramEsts(8),...
        paramEsts(9))));
end