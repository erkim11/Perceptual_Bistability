% initialize variables
s1 = zeros(1, 40001);
s2 = zeros(1, 40001);
s1(1) = 0.51;
s2(1) = 0.49;

% initialize adaptation and noise
a1 = zeros(1, 40001);
a2 = zeros(1, 40001);
n1 = zeros(1, 40001);
n2 = zeros(1, 40001);

% parameters
dt = 0.1; % ms
g = 20; % adaptation effect
sigma = 0; % zero noise
i = 1; % index

for t = 0:dt:4000
    % update synaptic action variables
    s1(i + 1) = s1(i) + dt * ds_dt(s1(i), s2(i), a1(i), n1(i), g);
    s2(i + 1) = s2(i) + dt * ds_dt(s2(i), s1(i), a2(i), n2(i), g);

    % update adaptation variables
    a1(i + 1) = a1(i) + dt * da_dt(a1(i), s2(i), n1(i), g);
    a2(i + 1) = a2(i) + dt * da_dt(a2(i), s1(i), n2(i), g);

    % update noise terms
    n1(i + 1) = n1(i) + dt * dn_dt(n1(i), sigma);
    n2(i + 1) = n2(i) + dt * dn_dt(n2(i), sigma);

    i = i + 1;
end

figure
figWidth = 1500;
figHeight = 400;
set(gcf, 'Position', [100, 100, figWidth, figHeight])
plot(s1)
hold on
plot(s2)
xlim([0 40000])
ylim([0 1])
xlabel('time step (0.1 ms)', 'FontSize', 16)
ylabel('firing rate', 'FontSize', 16)
title('g = 20', 'FontSize', 16)
legend('s1', 's2', 'FontSize', 16)
hold off

% no adaptation

% initialize variables
s1 = zeros(1, 40001);
s2 = zeros(1, 40001);
s1(1) = 0.5;
s2(1) = 0.5;

% initialize adaptation and noise
a1 = zeros(1, 40001);
a2 = zeros(1, 40001);
n1 = zeros(1, 40001);
n2 = zeros(1, 40001);

% parameters
dt = 0.1; % ms
g = 0; % adaptation effect
b = 7; % feedforward input
sigma = 20; % noise
i = 1; % index

for t = 0:dt:4000
    % update synaptic action variables
    s1(i + 1) = s1(i) + dt * ds_dt_b(s1(i), s2(i), a1(i), n1(i), g, b);
    s2(i + 1) = s2(i) + dt * ds_dt_b(s2(i), s1(i), a2(i), n2(i), g, b);

    % update adaptation variables
    a1(i + 1) = a1(i) + dt * da_dt(a1(i), s2(i), n1(i), g);
    a2(i + 1) = a2(i) + dt * da_dt(a2(i), s1(i), n2(i), g);

    % update noise terms
    n1(i + 1) = n1(i) + dt * dn_dt(n1(i), sigma);
    n2(i + 1) = n2(i) + dt * dn_dt(n2(i), sigma);

    i = i + 1;
end

figure
figWidth = 1500;
figHeight = 400;
set(gcf, 'Position', [100, 100, figWidth, figHeight])
plot(s1)
hold on
plot(s2)
xlim([0 40000])
ylim([0 1])
xlabel('time step (0.1 ms)', 'FontSize', 16)
ylabel('firing rate', 'FontSize', 16)
title('sigma = 20', 'FontSize', 16)
legend('s1', 's2', 'FontSize', 16)
hold off

% compute average time interval
mean_isis = zeros(1, 41);

j = 1;
for sigma = 0:0.5:20
    % initialize variables
    s1 = zeros(1, 60001);
    s2 = zeros(1, 60001);
    s1(1) = 0.51;
    s2(1) = 0.49;
    
    % initialize adaptation and noise
    a1 = zeros(1, 60001);
    a2 = zeros(1, 60001);
    n1 = zeros(1, 60001);
    n2 = zeros(1, 60001);
    
    % parameters
    dt = 0.1; % ms
    g = 20; % adaptation effect
    b = 7; % feedforward input
    i = 1; % index
    
    for t = 0:dt:6000
        % update synaptic action variables
        s1(i + 1) = s1(i) + dt * ds_dt_b(s1(i), s2(i), a1(i), n1(i), g, b);
        s2(i + 1) = s2(i) + dt * ds_dt_b(s2(i), s1(i), a2(i), n2(i), g, b);
    
        % update adaptation variables
        a1(i + 1) = a1(i) + dt * da_dt(a1(i), s2(i), n1(i), g);
        a2(i + 1) = a2(i) + dt * da_dt(a2(i), s1(i), n2(i), g);
    
        % update noise terms
        n1(i + 1) = n1(i) + dt * dn_dt(n1(i), sigma);
        n2(i + 1) = n2(i) + dt * dn_dt(n2(i), sigma);
    
        i = i + 1;
    end
    
    % isi distribution
    s_diff = s1 - s2;
    
    signs = sign(s_diff);
    diff_signs = diff(signs);
    switches = find(diff_signs);
    
    isi = diff(switches);
    mean_isis(j) = mean(isi);
    j = j + 1;
end

sigma = 0:0.5:20;

figure
figWidth = 900;
figHeight = 600;
set(gcf, 'Position', [100, 100, figWidth, figHeight])
plot(sigma, mean_isis)
xlabel('sigma', 'FontSize', 16)
ylabel('average time interval between switches (0.1 ms)', 'FontSize', 16)
hold off

% plot isi distribution over multiple sigma
figure
figWidth = 1500;
figHeight = 300;
set(gcf, 'Position', [100, 100, figWidth, figHeight])

subplot(1, 4, 1)
histogram(isi2, 20)
% ylim([0 56]);
xlabel('time intervals (0.1 ms)', 'FontSize', 16);
ylabel('count', 'FontSize', 16);
title('sigma = 2', 'FontSize', 16);

subplot(1, 4, 2)
histogram(isi4, 20)
% ylim([0 56]);
xlabel('time intervals (0.1 ms)', 'FontSize', 16);
ylabel('count', 'FontSize', 16);
title('sigma = 4', 'FontSize', 16);

subplot(1, 4, 3)
histogram(isi6, 20)
% ylim([0 56]);
xlabel('time intervals (0.1 ms)', 'FontSize', 16);
ylabel('count', 'FontSize', 16);
title('sigma = 6', 'FontSize', 16);
hold off

sigma = 0:0.5:20;
subplot(1, 4, 4)
plot(sigma, mean_isis)
xlabel('sigma', 'FontSize', 16)
ylabel('mean time intervals (0.1 ms)', 'FontSize', 16)
hold off

% compute STA
% no adaptation

% initialize variables
s1 = zeros(1, 400001);
s2 = zeros(1, 400001);
s1(1) = 0.5;
s2(1) = 0.5;

% initialize adaptation and noise
a1 = zeros(1, 400001);
a2 = zeros(1, 400001);
n1 = zeros(1, 400001);
n2 = zeros(1, 400001);

% parameters
dt = 0.1; % ms
g = 0; % adaptation effect
b = 7; % feedforward input
sigma = 20; % noise
i = 1; % index

for t = 0:dt:40000
    % update synaptic action variables
    s1(i + 1) = s1(i) + dt * ds_dt_b(s1(i), s2(i), a1(i), n1(i), g, b);
    s2(i + 1) = s2(i) + dt * ds_dt_b(s2(i), s1(i), a2(i), n2(i), g, b);

    % update adaptation variables
    a1(i + 1) = a1(i) + dt * da_dt(a1(i), s2(i), n1(i), g);
    a2(i + 1) = a2(i) + dt * da_dt(a2(i), s1(i), n2(i), g);

    % update noise terms
    n1(i + 1) = n1(i) + dt * dn_dt(n1(i), sigma);
    n2(i + 1) = n2(i) + dt * dn_dt(n2(i), sigma);

    i = i + 1;
end

% isi distribution
s_diff = s1 - s2;

signs = sign(s_diff);
diff_signs = diff(signs);

positive_elements = diff_signs > 0;
negative_elements = diff_signs < 0;

s1_sup_dom = zeros(size(diff_signs));
s1_sup_dom(positive_elements) = 1;
s2_sup_dom = zeros(size(diff_signs));
s2_sup_dom(negative_elements) = 1;
s1_dom_sup = zeros(size(diff_signs));
s1_dom_sup(negative_elements) = 1;
s2_dom_sup = zeros(size(diff_signs));
s2_dom_sup(positive_elements) = 1;

% compute sta for s1
[c_s1, lag_s1] = xcorr(n1, s1_dom_sup);
% c_s1_f = fliplr(c_s1);
mid = floor(length(lag_s1)/2);
final_c_s1 = c_s1(mid - 12000: mid + 12000) / sum(s1_dom_sup)
final_lag_s1 = lag_s1(mid - 12000: mid + 12000)

% compute sta for s2
[c_s2, lag_s2] = xcorr(n2, s2_sup_dom);
% c_s2_f = fliplr(c_s2);
mid = floor(length(lag_s2)/2);
final_c_s2 = c_s2(mid - 12000: mid + 12000) / sum(s2_sup_dom)
final_lag_s2 = lag_s2(mid - 12000: mid + 12000)

figure
figWidth = 1500;
figHeight = 500;
set(gcf, 'Position', [100, 100, figWidth, figHeight])

subplot(1, 2, 1)
plot(b * (1 + final_lag_s1), final_c_s1)
hold on
plot(b * (1 + final_lag_s2), final_c_s2)
xlim([-1000 1000])
xlabel('time to switch (0.1 ms)', 'FontSize', 16)
ylabel('stimulus', 'FontSize', 16)
legend('s1', 's2', 'FontSize', 16)
title('in-to-out direction', 'FontSize', 16)
hold off

% compute sta for s1
[c_s1, lag_s1] = xcorr(n1, s1_sup_dom);
% c_s1_f = fliplr(c_s1);
mid = floor(length(lag_s1)/2);
final_c_s1 = c_s1(mid - 12000: mid + 12000) / sum(s1_sup_dom)
final_lag_s1 = lag_s1(mid - 12000: mid + 12000)

% compute sta for s2
[c_s2, lag_s2] = xcorr(n2, s2_dom_sup);
% c_s2_f = fliplr(c_s2);
mid = floor(length(lag_s2)/2);
final_c_s2 = c_s2(mid - 12000: mid + 12000) / sum(s2_dom_sup)
final_lag_s2 = lag_s2(mid - 12000: mid + 12000)

subplot(1, 2, 2)
plot(b * (1 + final_lag_s1), final_c_s1)
hold on
plot(b * (1 + final_lag_s2), final_c_s2)
xlim([-1000 1000])
xlabel('time to switch (0.1 ms)', 'FontSize', 16)
ylabel('stimulus', 'FontSize', 16)
legend('s1', 's2', 'FontSize', 16)
title('out-to-in direction', 'FontSize', 16)
hold off
