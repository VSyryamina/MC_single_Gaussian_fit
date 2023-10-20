clear 
close all
%% The script for fitting the normalized dipolar trace by a single Gaussian distribution
%% Returns the best fit as (r_0, dr_0) pair and all (r_0, dr_0, loss) as table data_out

% example is given for trace #17
% where the ground truth is Gaussian distribution with r_0 = 3.96, dr_0 = 0.396:
% ground_truth = 1 / dr_0 /(2*pi)^0.5 * exp(-(r-r_0).^2/2/ dr_0^2);

%% load experimental trace
% the trace should be background-corrected and normalized to the modulation
% depth as described by von Hagens (2013) or Milov (2013)
% in the example data are the table with time t and normalized dipolar decay V
data = readtable('trace_17.txt');
t = data.t; V_exp = data.V;

%% define parameters
r = 1.8:0.01:6; % distance range, nm
ground_truth = 1 / (0.1*3.96)/(2*pi)^0.5 * exp(-(r-3.96).^2/2/ (0.1*3.96)^2);
pB = 1;  % excitation probability, also is known as lambda; pB == 1 for normalized dipolar trace
gAB = 1; % gAB == 1 if gA, gB are for free electrons; either one need to use gAB = gA*gB/g_free_electron^2

%% defining of the fitting parameters range in Gaussian parametric model
r_0_range = [3.6 4.20];
dr_0_over_r_0_range = [0.03 0.2];

[r_0_min, r_0_diap] = fit_range (r_0_range);
[dr_0_over_r_0_min, dr_ver_r_diap] = fit_range (dr_0_over_r_0_range);

%%
MC_trials = 5000; % number of trials
data_out.r_0 = zeros(1, MC_trials); data_out.dr_0 = zeros(1, MC_trials); data_out.loss = zeros(1, MC_trials);

parfor i = 1:MC_trials
    %% Create probe distance distribution function
    r_0_test = r_0_min+rand(1)*r_0_diap;   dr_0_test = r_0_test*(dr_0_over_r_0_min+rand(1)*dr_ver_r_diap);
    f_r_test = 1 / dr_0_test/(2*pi)^0.5 * exp(-(r-r_0_test).^2/2/dr_0_test^2);

    %% create dipolar decay
    % use DeerLab function dipolarsignal (https://jeschkelab.github.io/DeerLab/) (faster) 
    V_sim = dipolarsignal(t,r,f_r_test,pB)
    % or custom script
    % V_sim = create_DEER_decay (r, f_r_test, t, pB, gAB, 0);

    %% calculate mse between fit and experiment
    loss = sqrt(mean((V_sim-V_exp).^2));

    % keep data in memory for error-surface and lookin the best (r_0, dr_0) pair
    data_out(i).r_0 = r_0_test;
    data_out(i).dr_0 = dr_0_test;
    data_out(i).loss = loss;
end

[data_out, data_best] = treat_data_out(data_out);
best_fit = [mean(data_best.r_0), mean(data_best.dr_0)]

f_r_best = 1 / best_fit(2)/(2*pi)^0.5 * exp(-(r-best_fit(1)).^2/2/best_fit(2)^2);
V_best = create_DEER_decay(r, f_r_best, t, pB, gAB, 0);

%%
figure(1)
scatter(t, V_exp, 8, 'filled'); hold on
plot(t, V_best, 'linewidth', 2); hold on
legend('experimental decay','best fit'); legend boxoff
xlabel ('\itt, \mus');  ylabel ('\itV(t)'); 

figure(2)
plot(r, ground_truth, 'k', 'linewidth', 2); hold on
plot(r, f_r_best, 'linewidth', 2); hold on
legend('ground truth','best fit'); legend boxoff
xlabel ('\itr, nm');  ylabel ('\itf(r), nm^{-1}'); 

size = 12;
figure(3)
scatter3(data_out.r_0, data_out.dr_0, data_out.loss, size, data_out.loss, 'filled'); hold on
xlabel ('\itr_0, nm');  ylabel ('\it\deltar_0, nm'); 
view(0,90)
colormap jet
colorbar