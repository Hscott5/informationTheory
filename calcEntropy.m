function [H, bins] = calcEntropy(dat, type, bins)

if nargin<2
    type ='data';
end
if nargin<3
    bins=[];
end
nSimulations = 300; % for data simulations

dat = dat(:);
if numel(unique(dat))==1
    H = NaN;
    bins = [];
    return
end

switch type
    case 'data'
        if ~any(dat)
            H = NaN;
            return
        else
            % if no bin input was given
            if isempty(bins)
                binSize = 3.49*std(dat)*numel(dat)^(-1/3); % Scott 1979
                bins = min(dat)-binSize:binSize:max(dat)+binSize;
            end
            n = histcounts(dat, bins);
            n(n==0) = [];
            if ~isempty(n)
                p = n./sum(n);
                H = -sum(p.*log2(p));
            else
                H = NaN;
            end
        end
    case 'poiss'
        if ~any(dat)
            H = NaN;
            return
        else
            lambda = mean(dat); 
           
            % nSimulations for that lambda to get a better estimate
            s = poissrnd(lambda, [nSimulations,numel(dat)]);
            if isempty(bins)
                binSize = 3.49*std(s(:))*numel(dat)^(-1/3); % Scott 1979
                bins = min(dat,[],'all')-binSize:binSize:max(dat,[],'all')+binSize;
            end
            sn = cellfun( @(x) histcounts(x,bins), mat2cell(s, ones(nSimulations,1), numel(dat)), 'UniformOutput',false);
%             sn = cellfun(@(x) x(x~=0), sn, 'UniformOutput', false);
            pn = cellfun(@(x) x./sum(x), sn, 'UniformOutput',false);
            H = mean(cellfun(@(x) -nansum(x.*log2(x)),pn));
        end
    case 'theory'
        lambda = mean(dat);
        if lambda < 100
            k = 1:150;
            S= exp(-lambda).*sum((lambda.^k./factorial(k)).*log(factorial(k)));
            H = lambda*(1- log(lambda)) + S;
        else
            H = 0.5*log(2*pi*exp(1)*lambda)- (1/(12*lambda))-1/((24*lambda)^2);
        end
end


