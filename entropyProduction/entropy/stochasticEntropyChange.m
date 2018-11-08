% function [deltaS, edges, bincenters, vss] = stochasticEntropyChange(tseries, dt, dbin, d)
%
% function to calculate the stochastic entropy change along a trajectory
% defined by the input array.
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
%
% Returns
% -------
% deltaS : float
%     total entropy production
% edges : cell array
%     cell array with position of edges of n-dimensional histogram that is created to
%     discretize space. edges{n} gives edges in n-th dimension.
% bincenters : cell array
%     cell array with positions of centers of bins. bincenters{n} gives bin centers in n-th dimension
% vss : array
%     numerical array with velocity vectors over discretized space. First D dimensions index which
%     bin the velocity belongs to, the D+1 dimension contains all the vector components at that point
%
function [deltaS, edges, bincenters, vss] = stochasticEntropyChange(tseries, dt, dbin, d)
    [T,D] = size(tseries);

    if nargin < 4
        d = 'cov';
    end

    if isa(d, 'string') || isa(d, 'char')
        dInv = inv(cov(diff(tseries)) ./ dt);
    elseif isa(d, 'cell')
        d = rodTransverseDiffusion(d{:});
        dInv = inv(d);
    elseif isa(d, 'numeric')
        dInv = inv(d);
    end
    % Get the edges to discretize space
    edges = cell(1,D);
    nbin = cell(1,D);
    bincenters = cell(1,D);
    for ii = 1:D
        nbin{ii} = ceil((max(tseries(:,ii)) - min(tseries(:,ii))) / dbin);
        edges{ii} = linspace(min(tseries(:,ii)), max(tseries(:,ii)), nbin{ii} + 1);
        centers = edges{ii} + dbin/2;
        bincenters{ii} = centers(1:end-1);
    end

    deltaS = zeros(T-2,1);
    vss = zeros(nbin{:}, D);
    visitedStates = zeros(nbin{:}, D);

    % get steady state velocity
    for jj = 2:T-1
        % Get current state
        c = cell(1, D);
        currentState = tseries(jj,:);
        [currentCounts,~,~,~] = histcountsn(currentState, edges{:});
        visitedStates = visitedStates + currentCounts;
        [c{:}] = ind2sub(size(currentCounts), find(currentCounts));

        % calculate velocity
        v = (tseries(jj+1, :) - tseries(jj-1, :)) ./ (2 * dt);
        % add to velocity field over discrete space need to reshape
        % v to be [1,1,1,1,..., D], where the ones are repeated D times
        vss(c{:}, :) = vss(c{:}, :) + reshape(v, [ones(1, D), D]);
    end
    % get average
    vss = vss ./ visitedStates; %(T-2);

    % Calculate entropy produced at each time point
    for jj = 2:T-1
        % Get current state
        c = cell(1, D);
        currentState = tseries(jj,:);
        [currentCounts,~,~,~] = histcountsn(currentState, edges{:});
        [c{:}] = ind2sub(size(currentCounts), find(currentCounts));

        v = (tseries(jj+1, :) - tseries(jj-1, :)) ./ (2 * dt);
        deltaS(jj) = dt .* v * dInv * squeeze(vss(c{:}, :));
    end
end
