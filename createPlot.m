function createPlot(InMatrix1, type, info, Data)

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(InMatrix1,'Parent',axes1, 'LineWidth',2);
switch type
    case 'states'
        switch Data.Model_choice
            case 'BS'
                if size(InMatrix1,2) == 4
                    set(plot1(1),'DisplayName','Volatility Estimate');
                    set(plot1(2),'DisplayName','Rate Estimate');
                    set(plot1(3),'DisplayName','Implied Volatility');
                    set(plot1(4),'DisplayName','Rate Actual');
                elseif size(InMatrix1,2) == 2
                    set(plot1(1),'DisplayName','Volatility Estimate');
                    set(plot1(2),'DisplayName','Implied Volatility');
                end
            case 'md_RBF'
                set(plot1(1),'DisplayName','Mean 1 Estimate');
                set(plot1(2),'DisplayName','Mean 2 Estimate');
                set(plot1(3),'DisplayName','Mean 1 Actual');
                set(plot1(4),'DisplayName','Mean 2 Actual');
        end
        % Create ylabel
        ylabel({'state vector'});
        
        % Create title
        title({'Evolution of the state vector', info});
    case 'observation'
        set(plot1(1),'DisplayName','Measurement Estimate');
        set(plot1(2),'DisplayName','Measurement Actual');
        
        
        % Create ylabel
        ylabel({'Measurement(Z)'});
        
        % Create title
        title({'Comparing True and Estimanted Observation', info});
        set(gca,'FontSize',12);
        
end
% Create xlabel
xlabel({'time'});
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',12);


