function [isOutside] = checkFilPosition(x,y,xlims,ylims)
    xlo = xlims(1);
    xhi = xlims(2);
    ylo = ylims(1);
    yhi = ylims(2);

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
