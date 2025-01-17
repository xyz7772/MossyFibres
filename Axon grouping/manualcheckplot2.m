function manualcheckplot2(dff0_r, groups_cc, groups_ldr)

    groupidx = 1;
    totalgroups = length(groups_ldr);
    groups = [];
    tickStatus = cell(size(groups_ldr)); 
    check = [];
    fig = [];
    historyStack = [];
    
    createFig();
    
    function createFig()
        if groupidx > totalgroups
            handleUngrouped();
            return;
        end

        if ishandle(fig)
            close(fig);
        end

        m = groups_cc{groupidx};
        if isempty(m) || length(m) == 1 
            groups{groupidx} = m;
            groupidx = groupidx + 1;
            createFig();
            return;
        elseif length(groups_ldr) >= groupidx && isequal(sort(m), sort(groups_ldr{groupidx}))
            groups{groupidx} = m;
            groupidx = groupidx + 1;
            createFig();
            return;
        end

        if ~isempty(historyStack) && historyStack(end) ~= groupidx
            historyStack(end+1) = groupidx;
        elseif isempty(historyStack)
            historyStack(end+1) = groupidx;
        end
        
        fig = figure('Visible', 'on', 'Position', [100, 50, 1200, 750], 'CloseRequestFcn', @closeFigure, 'WindowStyle', 'normal', 'WindowState', 'maximized');
        drawnow;
        
        check = [];
        for plot_i = 1:length(m)
            ax = subplot(length(m), 1, plot_i);
            plotColor = 'k';
            if ~isempty(groups_ldr) && length(groups_ldr) >= groupidx && ~ismember(m(plot_i), groups_ldr{groupidx})
                plotColor = 'r';
            end
            
            plot(ax, dff0_r(m(plot_i), :), plotColor);
            set(ax, 'LineWidth', 0.5, 'FontSize', 6, 'FontName', 'Helvetica', 'Box', 'off', 'TickDir', 'out', 'Ticklength', [0.002, 0.002], 'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off', 'Xtick', []);
            grid(ax, 'minor');
   
            pos = get(ax, 'Position');
            boxWidth = 0.04;
            boxHeight = 0.05;
            boxLeft = pos(1) + pos(3) - boxWidth;
            boxBottom = pos(2) + pos(4)/2 - boxHeight/2;
            check(plot_i) = uicontrol('Style', 'checkbox', 'String', ['#' num2str(m(plot_i))], 'Units', 'normalized', 'Position', [boxLeft, boxBottom, boxWidth, boxHeight], 'Value', plot_i == 1 || plotColor == 'k');
        end

        assignin('base', 'manual_groups', groups);
        disp(groups)

        uicontrol('Style', 'text', 'String', sprintf('Group: %d of %d', groupidx, totalgroups), 'Position', [0.1 0.1 180 40], 'FontSize', 15);
        uicontrol('Style', 'pushbutton', 'String', 'Next', 'Position', [350 20 100 40], 'Callback', {@nextBtn});
        uicontrol('Style', 'pushbutton', 'String', 'Back', 'Position', [500 20 100 40], 'Callback', {@backBtn});
        uicontrol('Style', 'pushbutton', 'String', 'Close All', 'Position', [650 20 100 40], 'Callback', @closeAllBtn);
    end

    function nextBtn(src, event)    
        tickStatus = arrayfun(@(x) get(x, 'Value'), check);
        selected = find(tickStatus == 1);
        groups{groupidx} = groups_cc{groupidx}(selected);
        close(gcf);
        groupidx = groupidx + 1;
        createFig();
    end

    function backBtn(src, event)
        if length(historyStack) > 1
            historyStack(end) = [];
            groupidx = historyStack(end);
            historyStack(end) = [];
            createFig();
        end
    end

    function closeFigure(src, event)
        delete(gcf);
    end

    function closeAllBtn(src, event)
        delete(findall(0, 'Type', 'figure'));
    end

    function handleUngrouped()
    disp('Final manual groups:');
    disp(groups);
    end

end
