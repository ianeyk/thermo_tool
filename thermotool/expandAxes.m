function expandAxes(scale) % must have current fig selected
    if ~exist('scale', 'var')
        scale = 1.1
    end
    xlim((xlim - sum(xlim ./ 2)) .* scale + sum(xlim ./ 2));
    ylim((ylim - sum(ylim ./ 2)) .* scale + sum(ylim ./ 2));
end