function [x_reduced, y_reduced] = reduce_to_height(x, y, lims)

% [x_reduced, y_reduced] = reduce_to_height(x, y, lims)
% 
% For seismograms only, all values of each trace are near trace number (F2PData2 function)
% 
% Geser Dugarov 2016

    % If the data is already small, there's no need to reduce.
    if max(y{end}(:))-min(y{1}(:)) <= lims(2)-lims(1)
        x_reduced = x;
        y_reduced = y;
        return;
    end

    % Reduce the data to the new Y axis size.
    bot_bound = max([1,ceil(lims(1))]);
    upp_bound = min([size(y, 2),floor(lims(2))]);
    x_reduced = cell(1, upp_bound-bot_bound+1);
    y_reduced = cell(1, upp_bound-bot_bound+1);
    for k = bot_bound:upp_bound
        x_reduced{k-bot_bound+1} = x{k}(:);
        y_reduced{k-bot_bound+1} = y{k}(:);
    end
end