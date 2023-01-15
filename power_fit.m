function [fitresult, gof] = power_fit(t1, t2)

[xData, yData] = prepareCurveData( t1, t2 );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxFunEvals = 6000;
opts.MaxIter = 4000;
opts.Robust = 'Bisquare';
opts.StartPoint = [937.568059384175 -0.308372447724913];
opts.TolFun = 1e-07;
opts.TolX = 1e-07;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 't2 vs. t1', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 't1', 'Interpreter', 'none' );
% ylabel( 't2', 'Interpreter', 'none' );
% grid on
end