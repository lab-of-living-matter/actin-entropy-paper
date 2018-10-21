% function [isOutside] = checkFilPosition(x,y,xlims,ylims)
%
% Function that checks whether a given set of position is
% beyond some bounding box. Use to check whether a cytosim filament
% approaches the periodic boundary
%
% Parameters
% ----------
% x, y: array
%     x and y positions to check
% xmax, ymax: float
%     maximum absolute value of x and y. For example, if simulation
%     box goes from [-10, 10] in x, and [-15,15] in y, xmax=10, ymax=15
% percent : float
%     percentage of xmax and ymax to allow filaments to traverse. Write as
%     decimal, so 95 percent is 0.95
%
% Returns
% -------
% isOutside: bool
%     boolean that says whether the x,y positions go outside
%     +/- (percent) * (xmax or ymax)
%
% Created by Ian Linsmeier, 06/07/2018
function [isOutside] = checkFilPosition(x, y, xmax, ymax, percent)
    xlo = -xmax .* percent;
    xhi = xmax .* percent;
    ylo = -ymax .* percent;
    yhi = ymax .* percent;

    isOutside = logical(0);
    xloTest = any(x < xlo);
    isOutside = xloTest + isOutside;
    if isOutside == 1
        return
    end

    xhiTest = any(x > xhi);
    isOutside = xhiTest + isOutside;
    if isOutside == 1
        return
    end

    yloTest = any(y < ylo);
    isOutside = yloTest + isOutside;
    if isOutside == 1
        return
    end

    yhiTest = any(y > yhi);
    isOutside = yhiTest + isOutside;
    if isOutside == 1
        return
    end
end
