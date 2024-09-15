% In this script, six values are closely approximated via golden section search
%  - α value for which correlation coefficient ρ is strongest between P_bp and T_{a}^1
%  - α value for which ρ is strongest between P_bp and T_{a}^2
%  - α value for which ρ is strongest between Δ_{h0f} and T_{a}^1
%  - α value for which ρ is strongest between Δ_{h0f} and T_{a}^2
% T_{a}^1 = First General Temperature Indices
% T_{a}^2 = Second General Temperature Indices
% Additionally, curves for ρ against α near these values are plotted in 4 figs

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console
pkg load statistics;

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2; % the line thickness in the graph
fontSize = 16; % the font size in the graph
saveToFile = true; % Set to true to auto-save plots
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
numData = size(expData, 2);       % two (ρ_{bp}, Δ_{h0f})
numIndices = size(getIndexFns,2); % two (T_{a}^1  & T_{a}^2)
numCurves = numData*numIndices;   % four (4 curves)

% All x in visible ranges (both plots - near and far)
xnear = [linspace(-0.1, 0.45, 800); linspace(-0.8, 0.5, 800)];
xfar = linspace(-20,20,800); % xfar is same for both ρ_{bp} and Δ_{h0f}

% Do the same procedure for each experimental data i.e., ρ_{bp} and Δ_{h0f}
for ii = 1:numData
  % This figure is for the zoomed-in plot
  figure(numData+ii); hold on;

  % Indicate inverval for which R_a is the better indicator (before corr-curves)

  % WARNING: these xmeet1-ymeet2 values are hardcoded, computed separately
  %          ... These are coordinates where ρ-α of {Δ_{h0f} or E}-T_{a}^1
  %          ... intersects ρ-α of {Δ_{h0f} or E}-T_{a}^2
  xmeet1 = [
    0.1349; % for ρ_{bp}
    -0.12142   % for Δ_{h0f}
  ](ii); % of these 2, either the first value or second is used, depending on ii
  ymeet1 = [
    0.99692; % for ρ_{bp}
    0.94592  % for Δ_{h0f}
  ](ii);
  xmeet2 = [
    0.00012;    % for ρ_{bp}
    -0.00010081  % for Δ_{h0f}
  ](ii);
  ymeet2 = [
    0.99576; % for ρ_{bp}
    0.94512  % for Δ_{h0f}
  ](ii);

  ybox = [
    0.9972; % for figure 3 (blue dotted line)
    0.9465  % for figure 4 (blue dotted line)
  ](ii);
  % Plot the blue dashed box (before drawing the curves so it appear beneath)
  plot([xmeet1 xmeet1 0 0], [0 ybox ybox 0], '--b', 'LineWidth', lineWidth);

  yend = 0.995; % <-- to be assigned some value later for adjusting visible range

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between ρ_{bp}/Δ_{h0f} (depending on ii) with specified α
    %                                and T_{a}^1/T_{a}^2 (depending on n)
    get_indices_vals = @(alpha) getIndexFns{n}(alpha)(!isnan(expData{1,ii}));
    ccFn = @(alpha) corrcoef(
      get_indices_vals(alpha), % Either T_{a}^1  or T_{a}^2
      expData{1,ii}(!isnan(expData{1,ii})) % ρ_{bp}/Δ_{h0f}
    )(1,2);

    % generate corresponding y values
    ynear = arrayfun(ccFn, xnear(ii,:));
    yfar = arrayfun(ccFn, xfar);

   % Compute peak values via. golden section search, and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakAlpha = mean(GoldenSectionSearch_Maximum(ccFn, xnear(1), xnear(end), 1e-15));
    peakCorrCoeff = ccFn(peakAlpha);

    % Display peak alpha and peak correlation coefficient in the command window
    disp(['Peak Alpha: ', num2str(peakAlpha)]);
    disp(['Peak Correlation Coefficient: ', num2str(peakCorrCoeff)]);

    % Generate curve label                  [ρ_{bp}/Δ_{h0f}]         [T_{a}^1/T_{a}^2 ]
    curveLabels{n} = sprintf("%s and %s_a", expData{2,ii}, indexName{n});

    figure(ii); % One zoomed-out plot for each expData
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;

    figure(numData+ii); % Each expData's set of curves has a zoomed-in plot
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curves(n) = plot(xnear(ii,:), ynear, '-', 'LineWidth', lineWidth);

    % Show the peak in the plot: draw indicator lines & display coordinates
    plot([peakAlpha peakAlpha xnear(ii,1)], [0 peakCorrCoeff peakCorrCoeff],
         '--k', 'LineWidth', lineWidth/2); % Black dashed indicator lines
    text(peakAlpha, peakCorrCoeff,
        {'', sprintf("(−%s, %s)", as_4_dp_str(abs(peakAlpha)), as_4_dp_str(peakCorrCoeff))},
        'VerticalAlignment', 'bottom');
        % Negative sign entered manually here to bypass the default
        % ... usage of "hypen-minus" instead of "minus" (− vs -)

    yend = max(yend, ynear(end)); % y value to be used as visible y lower bound
  end

  % Mark and write on the plot the limits of alpha where T_{a}^1 is better than T_{a}^2
  plot([xmeet1 xmeet2], [ymeet1 ymeet2], '*b',
       'MarkerSize', 16, 'LineWidth', lineWidth/1.5); % Mark with blue asterisks
  text(xmeet1, ymeet1, {'', sprintf(" (−%s, %s)", as_4_dp_str(abs(xmeet1)), as_4_dp_str(ymeet1))},
       'VerticalAlignment', 'top', 'Color', [0, 0, 0.8]); % Write blue text
  text(xmeet2, ymeet2, {'', sprintf("(0, %s) ", as_4_dp_str(ymeet2))},
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', [0, 0, 0.8]);

  % Label this expData's zoomed-in plot
  xlabel('α');
  ylabel('ρ');
  leg = legend(curves, curveLabels); % curves contains all drawn "xnear" curves
  set(leg, 'location', "southwest"); % the location of the legend box

  % Ensure yend is less than ybox + ybox_space before setting axis limits
  ybox_space = [
      0.000005; % spacing between upper y-value with the blue dotted lines (figure 3)
      0.0005  % spacing between upper y-value with the blue dotted lines (figure 4)
  ](ii);
  y_upper_limit = ybox + ybox_space;

  % Check if yend is less than y_upper_limit, otherwise, adjust it
  if yend >= y_upper_limit
      yend = y_upper_limit - 0.01; % Ensure yend is slightly less than the upper limit
  end

  axis([xnear(ii,1) xnear(ii,end) yend y_upper_limit]); % Enforce figure's visible range
  drawnow;

  % Label the zoomed-out plot
  figure(ii);
  xlabel('α');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels); % curvesFar contains all drawn "xfar" curves
  set(leg, 'location', "southeast");

  if ii==2
    set(leg, 'location', "southeast");

  end

  hold off;
end

for ii = 1:4
  % Replace hyphens with minuses on negative axes
  figure(ii);
  xticklabels(strrep(xticklabels,'-','−'));
  yticklabels(strrep(yticklabels,'-','−'));

  % Set all fontsizes to size specified early in the script
  set(findall(figure(ii),'-property','FontSize'),'FontSize', fontSize)
end

if saveToFile
  saveas(figure(1), "01_comb_ccurves_bp_indices_FAR.png");
  saveas(figure(2), "01_comb_ccurves_DH_indices_FAR.png");
  saveas(figure(3), "01_comb_ccurves_bp_indices_NEAR.png");
  saveas(figure(4), "01_comb_ccurves_DH_indices_NEAR.png");
end

