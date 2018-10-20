% function [logPRatio, m] = dft(tSeries, dt, dbin, D, dsMax, nbins)
%
% Function to test whether a time series obeys the detailed fluctuation
% theorem (dft). It does if log(p(+S)/p(-S)) =~ S, wher ep(s) is the
% probability of observing a time series produce some amount of entropy, s
%
% Parameters
% ----------
% tseries : array
%     NxD matrix, each column is a time series in one variable
% dt : float
%     time interval between each row of data in tseries
% dbin : float
%     Size of bins (in tseries units), equal in all directions
% d : string, float, or cell
%     Diffusion constant. If d = 'cov', d = cov(diff(tseries)) ./ dt. If d is a cell,
%     d{:} are the inputs to rodTransverseDiffusion. Else can be specified as a constant
% dsMax : float
%     maximum entropy change used to test dft. Can only test a small region around zero
%     because destroying large amounts of entropy is exponentially difficult
% nbins : int
%     Number of bins to split p(s > -dsMax & s < dsMax) into
%
% Returns
% -------
% logPRatio : array
%     (nbins+1) array of log(p(+S)/p(-S))
% m : float
%     slope of linear fit to logPRatio vs. s
function [logPRatio, m] = dft(tSeries, dt, dbin, D, dsMax, nbins)

    if ~mod(nbins,2)
        error('nbins needs to be odd')
    end

    deltaS = stochasticEntropyChange(tSeries, dt, dbin, D);

    if isempty(find(deltaS >= dsMax)) || isempty(find(deltaS <= -dsMax))
        warning('dsMax given exceeds range of deltaS. Pick a smaller dsMax');
        logPRatio = NaN;
        m = NaN;
        return
    end

    [N, edges] = histcounts(deltaS(deltaS >= -dsMax & deltaS <= dsMax), linspace(-dsMax, dsMax, nbins + 1));

    if nnz(N) ~= numel(N)
        error('Got empty bins in histogram, make bins larger by choosing a smaller number for nbins')
    end

    logPRatio = log(N(1:numel(edges)/2) ./ flip(N(numel(edges)/2:end)));
    dftRangeDS = linspace(0, dsMax, numel(logPRatio));

    linearCoeffs = polyfit(dftRangeDS(~isinf(logPRatio)), logPRatio(~isinf(logPRatio)), 1);
    m = linearCoeffs(1);
end
