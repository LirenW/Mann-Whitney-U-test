function [ P ] = MWtest( A, B, flag, N, NLAT, NLON, NTIME )
%
%   (Bootstrap) Mann Whitney U test (also called the
%       Mann Whitney Wilcoxon (MWW), Wilcoxon rank-sum test,
%       or Wilcoxon Mann Whitney test)
%   ATTENTION: Parallel Computing Toolbox is NEEDED
%       Further examples and solutions in
%       https://ww2.mathworks.cn/help/distcomp/index.html
%
%     Inputs:
%      A, B - Input data, should be same dims with [lon lat time]
%      flag - Use bootstrap or not. True for open the bootstrap
%               resampling process.
%         N - The pairs of bootstrap resampling.
%     Outputs:
%       out - p-value of Mann Whitney U test. For two-tail
%               MW-U test, p < 0.025 and p > 0.075 means 
%               that pass the MW-U 95% significant level test.
%               -999 is the filling value.
%
%   Example:
%
%   out = MWtest(A,B,true,1000);

%  Created by Liren
%  Aug 2018
%
%  References:
%
% Efron, B. & Tibshirani, R. An Introduction to the Bootstrap. (Chapman and Hall,
%                 1993).
% Koirala, S., Hirabayashi, Y., Mahendran, R. and Kanae, S.: Global assessment of 
%                 agreement among streamflow projections using CMIP5 model outputs, 
%                 Environ. Res. Lett., 9(6), 064017, doi:10.1088/1748-9326/9/6/064017, 
%                 2014.

    P       = ones(NLON, NLAT) * -999;
    if flag == true
        parfor ix = 1 : NLON
            for iy = 1 : NLAT
                if isnan(A(ix, iy, :)) == 1 | isnan(B(ix, iy, :)) == 1
                    continue
                end

                X             = reshape(A (ix, iy,: ), NTIME, 1);
                Y             = reshape(B (ix, iy,: ), NTIME, 1);
                muX           = nanmean(X);
                muY           = nanmean(Y);
                muZ           = nanmean([X;Y]);

                [~, ~, stats] = ranksum(X, Y);
                U0            = stats.ranksum;

                Xt            = bootstrap(X, N);
                Yt            = bootstrap(Y, N);
                X2            = Xt + muZ - muX;
                Y2            = Yt + muZ - muY;

                for j = 1 : N
                    [~, ~, stats] = ranksum(X2(:, j), Y2(:, j));
                    Uj(j)         = stats.ranksum;
                end

                Uj      = sort(Uj);
                top     = find(Uj <= U0, 1, 'last' );
                bot     = find(Uj >= U0, 1, 'first');

                R0        = (top + bot)/2;
                P(ix, iy) = (R0 - 0.4) / (N + 0.2);
                [ix iy P(ix, iy)]

            end
        end
    else
        parfor ix = 1 : NLON
            for iy = 1 : NLAT
                if isnan(A(ix, iy, :)) == 1 | isnan(B(ix, iy, :)) == 1
                    continue
                end

                X             = reshape(A (ix, iy,: ), NTIME, 1);
                Y             = reshape(B (ix, iy,: ), NTIME, 1);
                muX           = nanmean(X);
                muY           = nanmean(Y);
                muZ           = nanmean([X;Y]);

                p = ranksum(X, Y);
            end
        end
    end
end
