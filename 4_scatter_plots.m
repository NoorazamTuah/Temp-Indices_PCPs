% In this script, six scatter plots are generated for the following pairs
%  - P_bp and T_{a}^1
%  - P_bp and T_{a}^2
%  - Δ_{h0f} and T_{a}^1
%  - Δ_{h0f} and T_{a}^2
% ... with their respective regression lines.

close all; % Before drawing, close any figures already opened
clear;     % Clear all variables

% CONFIG: Line width and font fize for each curve in drawn figures
lineWidth = 2;
fontSize = 26;
% Save plots to images? Set to true if yes
saveToFile = false;
% Note: Dimensions of resulting images are scaled according to each window size.
%       To control the dimensions, after running the whole script, resize each
%       ... figure window and then run only the saveas functions
%       ... manually, by selection, at the end of this script

% Utility: Below is a function to round-off to 4 decimal places | returns string
%          Need to use this function as round(X,4,Type) does not exist in Octave
%          ... and sprintf("%.04f",X) does not round properly for some numbers.
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Cell containing Entropy and Heat Capacity of 22 lower benzenoids
expData = {reshape([% Boiling point
   80.1 218 338 340 431 425 429 440 496 493 497
   547  542 535 535 531 519 590 596 594 595 393
   ]', 22, 1), reshape([ % Enthalpy of Formation
   75.2 141   202.7 222.6 271.1 277.1 275.1 310.5 296 289.9 319.2
   323  301.2 348   335   336.3 336.9 296.7 375.6 366 393.3 221.3
]', 22, 1),
   "ρ_{bp}", "Δ_{h0f}" % Their labels
};

% 22 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6 6 7 6 8 7  9 6  7  8 8 6  7  9  8  8  9  6  8  8  9  6
  0 4 6 8 8 10 6 12 10 8 8 12 10 10 12 12 10 12 12 12 10 8
  0 1 3 2 5 4  6 3  7  8 8 9  10  7 6  6  7  12  9  9 10 5
]'; % Used for computing indices based on edge endpoint degree partitions

% Define nu_exp here before using it in getIndexFns Lower 22 Bhs (nu_exp: number of vertices)
nu_exp = [
    6 10 14 14 18 18 18 18 20 20 20 22 22 22 22 22 22 24 24 24 24 16
]';

% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(alpha, nu_exp) | n=1:T_{a}^1 , n=2:T_{a}^2, a = alpha
getIndexFns = {
    @(alpha) sum(d_f .* [(4 ./ (nu_exp - 2)), ((5 .* nu_exp - 12) ./ ((nu_exp - 2) .* (nu_exp - 3))), (6 ./ (nu_exp - 3))].^alpha, 2); % First General Temperature Indices (T_{a}^1)
    @(alpha) sum(d_f .* [(4 ./ ((nu_exp - 2).^2)), (6 ./ ((nu_exp - 2) .* (nu_exp - 3))), (9 ./ ((nu_exp - 3).^2))].^alpha, 2); % Second General Temperature Indices (T_{a}^2)
}';
% Cell containing their respective labels
indexName = {"T^1", "T^2"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % two
numIndices = size(getIndexFns,2); % two
numCurves = numData*numIndices;   % four

for edn = 1:numData % edn = experimental data number | 1=bp, 2=ΔH
  for fnn = 1:numIndices % fnn = function number | 1=T^1, 2=T^2
    ccFn = @(alpha) corrcoef( % Gets corrcoef between benzenoid prop and index
      getIndexFns{fnn}(alpha)(!isnan(expData{1,edn})),
      expData{1,edn}(!isnan(expData{1,edn}))
    )(1,2);

    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -5, 5, 1e-15));

    % Regression line i.e., y = mx + b;
    model = polyfit(getIndexFns{fnn}(peakAlpha), expData{1,edn},1);
    m = model(1); b = model(2);         % For the regression line
    x = [getIndexFns{fnn}(peakAlpha)(1) max(getIndexFns{fnn}(peakAlpha))];
    y = m*x + b;

    % Scatter plot
    this_figure = figure(3*(edn-1)+fnn); hold on;
    regLine = plot(x, y, '-', 'LineWidth', lineWidth);
    points = plot(getIndexFns{fnn}(peakAlpha), expData{1,edn}, '*', 'MarkerSize', 8, 'LineWidth', lineWidth/2);
    bestIndexLabel = sprintf("%s_{−%s}", indexName{fnn}, as_4_dp_str(abs(peakAlpha)));
    pointsLabel = sprintf("%s and %s", bestIndexLabel, expData{2,edn});

    % To find the standard error of fit
    % Calculate residuals
    residuals = expData{1, edn} - (m * getIndexFns{fnn}(peakAlpha) + b);
    % Calculate the sum of squared residuals (SSE)
    sse = sum(residuals.^2);
    % Calculate degrees of freedom (df)
    n = length(expData{1, edn});
    df = n - 2;  % Two parameters estimated: slope and intercept
    % Calculate the standard error of fit (SE)
    se = sqrt(sse / df);
    fprintf('Standard Error of Fit for %s and %s: %.4f\n', ...
    bestIndexLabel, expData{2, edn}, se);


    % Label the scatter plot
    title(sprintf('between %s and %s', expData{2,edn}, bestIndexLabel));
    xlabel(bestIndexLabel);
    ylabel(sprintf('%s', expData{2,edn}));
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    leg = legend("Regression Line", "Actual Data");
    set(leg, 'location', "southeast");

    % Change the font size to size set in the start of the script
    set(findall(this_figure,'-property','FontSize'),'FontSize', fontSize)

    drawnow;
  end
end

if saveToFile
  % Also enforce more suitable axes for figures 2 and 5
  saveas(figure(1), "03_scatter_E_T^1.png");
  figure(2);
  axis([0.035 0.1125 250 600]);
  saveas(figure(2), "03_scatter_E_T^2.png");
  saveas(figure(3), "03_scatter_DH_T^1.png");
  figure(4);
  axis([0.2 0.9625 50 400]);
  saveas(figure(4), "03_scatter_DH_T^2.png");
end
