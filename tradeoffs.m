%% tradeoffs.m
% This code plots the trade-offs between social distancing, testing, and
% contact tracing as defined by the analytical expression in the manuscript
% "Simple control for complex pandemics.
%
%% Create plot of trade-offs

% Establish parameters
clear all
close all

d = 10;

v = 0:0.01:1.5;

c1 = 0;
R0res1 = (1 + v*(d-1))./(1 - c1*v);

c2 = 0.5;
R0res2 = (1 + v*(d-1))./(1 - c2*v);

c3 = 0.8;
R0res3 = (1 + v*(d-1))./(1 - c3*v);

% Plot results

figure(1)
box on
hold on

area(v, R0res3, -1, 'FaceColor', 'k', 'FaceAlpha', 0.1)
area(v, R0res2, -1, 'FaceColor', 'k', 'FaceAlpha', 0.2)
area(v, R0res1, -1, 'FaceColor', 'k', 'FaceAlpha', 0.3)
ylabel('$$\mathcal{R}_0 =$$ basic reproduction number', 'Interpreter', 'latex', 'Color', 'k')
xlabel('$$\nu =$$ daily testing fraction', 'Interpreter', 'latex', 'Color', 'k')
set(gca, 'FontSize', 20)
set(gca, 'TickLabelInterpreter', 'latex')

% Label regimes and literature data
text(0.75, 6, 'Controlled', 'Color', 'k', 'FontSize', 20, 'Interpreter', 'latex') 
text(0.1, 8, 'Uncontrolled', 'Color', 'k', 'FontSize', 20, 'Interpreter', 'latex') 

delta = 0.4;
text(0.02, 2.2+delta, 'Wuhan', 'Color', 'k', 'FontSize', 14, 'Interpreter', 'latex') 
text(0.02, 4.5+delta, 'Global', 'Color', 'k', 'FontSize', 14, 'Interpreter', 'latex') 
text(0.02, 11+delta, 'Diamond Princess', 'Color', 'k', 'FontSize', 14, 'Interpreter', 'latex') 

ylim([0 15])
xlim([0 1])

plot([-1 1.5], [2.2 2.2], '-k', 'LineWidth', 1.5)
plot([-1 1.5], [4.5 4.5], '--k', 'LineWidth', 1.5)
plot([-1 1.5], [11 11], ':k', 'LineWidth', 1.5)

set(gca, 'Layer', 'top')