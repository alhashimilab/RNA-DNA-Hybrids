clear;
clc;

%% loading all parameters

parent_path = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-Figures-Main-and-SI';
filename = '2023-NMR-RD-params.xlsx';
filepath = sprintf('%s/%s',parent_path,filename);

params_titles = {'pGS_pH','pT_pH','pA_pH',...
    'kES1_pH','kES2_pH','kMinor_pH',...
    'kGStES1_pH','kES1tGS_pH','kGStES2_pH','kES2tGS_pH','kES1tES2_pH','kES2tES1_pH',...
    'pGS_pH_err','pT_pH_err','pA_pH_err',...
    'kES1_pH_err','kES2_pH_err','kMinor_pH_err',...
    'kGStES1_pH_err','kES1tGS_pH_err','kGStES2_pH_err','kES2tGS_pH_err','kES1tES2_pH_err','kES2tES1_pH_err'};

%% Exchange parameters

% Use this script to load raw data from the parent excel file, extract
% populations and kinetic rates and plot bar graphs for Hybrids (main) and
% DNA and RNA with CGC sequence context (SI)

%% RNA-DNA Hybrid NMR data

% path to fitted R1rho data dGrU, pH 7.4: /Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-RD-Figures/2023-RD-RawData/dGrU_pH7p4/global/3state-triangle-shared
% path to fitted R1rho data dTrG, pH 7.4: /Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-RD-Figures/2023-RD-RawData/dTrG_pH7p4/global/2state-shared
% path to fitted R1rho data dGrU, pH 7.8: /Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-RD-Figures/2023-RD-RawData/dGrU_pH8p0/global/3state-triangle-shared
% path to fitted R1rho data rGdT, pH 7.8: /Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-RD-Figures/2023-RD-RawData/dTrG_pH8p0/global/3state-triangle-shared

params_Hybrid_NMR_tbl   = readtable(filepath,'Sheet','Hybrids_NMR');

construct_name_Hybrid_all = params_Hybrid_NMR_tbl.Construct;
temperatures_Hybrid_all   = params_Hybrid_NMR_tbl.Temperature_K_;
pH_Hybrid_25C             = params_Hybrid_NMR_tbl.pH;

%% Figure 2E - Hybrid Exchange parameters

error_line_width = 1;
color_8p0 = [46/255,49/255,146/255];
y_ticks = 0:0.1:0.5;

params_Hybrid_25C         = params_Hybrid_NMR_tbl(:,[25:36,48:59]);
params_Hybrid_25C_to_save = [params_titles;table2cell(params_Hybrid_25C)];

Construct = construct_name_Hybrid_all;
pH        = pH_Hybrid_25C;

Hybrid_25C_tbl = cell2table(params_Hybrid_25C_to_save(2:end,:));
Hybrid_25C_tbl.Properties.VariableNames = params_Hybrid_25C_to_save(1,:);
Hybrid_25C_tbl = addvars(Hybrid_25C_tbl,pH,'Before','pGS_pH');
Hybrid_25C_tbl = addvars(Hybrid_25C_tbl,Construct,'Before','pH');

Hybrid_8p0_25C_tbl = Hybrid_25C_tbl([2,4],:);

%% Figure 2E - Hybrid Exchange parameters, pH 7.8, populations

pT_vals_8p0 = [Hybrid_8p0_25C_tbl.pT_pH(2,1);Hybrid_8p0_25C_tbl.pT_pH(1,1)];
pT_errs_8p0 = [Hybrid_8p0_25C_tbl.pT_pH_err(2,1);Hybrid_8p0_25C_tbl.pT_pH_err(1,1)];

pA_vals_8p0 = [Hybrid_8p0_25C_tbl.pA_pH(2,1);Hybrid_8p0_25C_tbl.pA_pH(1,1)];
pA_errs_8p0 = [Hybrid_8p0_25C_tbl.pA_pH_err(1,1);Hybrid_8p0_25C_tbl.pA_pH_err(1,1)];

bar_vals = [pT_vals_8p0.*100,pA_vals_8p0.*100];
bar_errs = [pT_errs_8p0.*100,pA_errs_8p0.*100];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p0;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p0;
hb(2).LineStyle = 'none';

axis([0.5 2.5 0.0 0.5]);

set(gca, 'FontSize', 24, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'Ytick',y_ticks,'TickLength',[0.015,0.015]);

ylabel('population [%]', 'FontSize', 36);

legend({'Tautomer', 'Anion'});

propedit;

%% Figure 2E - Hybrid Exchange parameters, pH 7.8, kex

ES1_vals_8p0 = [Hybrid_8p0_25C_tbl.kES1_pH(2,1);Hybrid_8p0_25C_tbl.kES1_pH(1,1)];
ES1_errs_8p0 = [Hybrid_8p0_25C_tbl.kES1_pH_err(2,1);Hybrid_8p0_25C_tbl.kES1_pH_err(1,1)];

ES2_vals_8p0 = [Hybrid_8p0_25C_tbl.kES2_pH(2,1);Hybrid_8p0_25C_tbl.kES2_pH(1,1)];
ES2_errs_8p0 = [Hybrid_8p0_25C_tbl.kES2_pH_err(2,1);Hybrid_8p0_25C_tbl.kES2_pH_err(1,1)];

bar_vals = [ES1_vals_8p0,ES2_vals_8p0];
bar_errs = [ES1_errs_8p0,ES2_errs_8p0];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p0;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p0;
hb(2).LineStyle = 'none';

axis([0.5 2.5 0.0 12000]);

set(gca, 'FontSize', 24, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'TickLength',[0.015,0.015]);

ylabel('k_{ex} [s^{-1}]', 'FontSize', 36);

%legend({'Tautomer', 'Anion'});

propedit;

%% Figure 2E - Hybrid Exchange parameters, pH 7.8, k_reverse

ES1_kb_vals_8p0 = [Hybrid_8p0_25C_tbl.kES1tGS_pH(2,1);Hybrid_8p0_25C_tbl.kES1tGS_pH(1,1)];
ES1_kb_errs_8p0 = [Hybrid_8p0_25C_tbl.kES1tGS_pH_err(2,1);Hybrid_8p0_25C_tbl.kES1tGS_pH_err(1,1)];

ES2_kb_vals_8p0 = [Hybrid_8p0_25C_tbl.kES2tGS_pH(2,1);Hybrid_8p0_25C_tbl.kES2tGS_pH(1,1)];
ES2_kb_errs_8p0 = [Hybrid_8p0_25C_tbl.kES2tGS_pH_err(2,1);Hybrid_8p0_25C_tbl.kES2tGS_pH_err(1,1)];

bar_vals = [ES1_kb_vals_8p0,ES2_kb_vals_8p0];
bar_errs = [ES1_kb_errs_8p0,ES2_kb_errs_8p0];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p0;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p0;
hb(2).LineStyle = 'none';

axis([0.5 2.5 0.0 12000]);

set(gca, 'FontSize', 24, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'TickLength',[0.015,0.015]);


ylabel('k_{reverse} [s^{-1}]', 'FontSize', 36);

%legend({'Tautomer', 'Anion'});

propedit;

%% Figure 2E - Hybrid Exchange parameters, pH 7.8, k_forward

ES1_kf_vals_8p0 = [Hybrid_8p0_25C_tbl.kGStES1_pH(2,1);Hybrid_8p0_25C_tbl.kGStES1_pH(1,1)];
ES1_kf_errs_8p0 = [Hybrid_8p0_25C_tbl.kGStES1_pH_err(2,1);Hybrid_8p0_25C_tbl.kGStES1_pH_err(1,1)];

ES2_kf_vals_8p0 = [Hybrid_8p0_25C_tbl.kGStES2_pH(2,1);Hybrid_8p0_25C_tbl.kGStES2_pH(1,1)];
ES2_kf_errs_8p0 = [Hybrid_8p0_25C_tbl.kGStES2_pH_err(2,1);Hybrid_8p0_25C_tbl.kGStES2_pH_err(1,1)];

bar_vals = [ES1_kf_vals_8p0,ES2_kf_vals_8p0];
bar_errs = [ES1_kf_errs_8p0,ES2_kf_errs_8p0];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p0;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p0;
hb(2).LineStyle = 'none';

axis([0.5 2.5 0.0 18.0]);

set(gca, 'FontSize', 24, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'TickLength',[0.015,0.015]);

ylabel('k_{forward} [s^{-1}]', 'FontSize', 36);

%legend({'Tautomer', 'Anion'});

propedit;

%% DNA and RNA NMR data

params_Nat2018_NMR_tbl   = readtable(filepath,'Sheet','DNA_RNA_previously_measured');

construct_name_Nat2018_all = params_Nat2018_NMR_tbl.Construct;
temperatures_Nat2018_all   = params_Nat2018_NMR_tbl.Temperature_K_;

%% Figure S5 - Exchange parameters, DNA CGC vs. RNA CGC

error_line_width = 1;
color_8p4 = [46/255,49/255,146/255];
y_ticks = 0:0.1:0.5;

indices_10C = find(ismember(temperatures_Nat2018_all,283.2));

params_Nat2018_10C         = params_Nat2018_NMR_tbl(indices_10C,[25:36,48:59]);
params_Nat2018_10C_to_save = [params_titles;table2cell(params_Nat2018_10C)];
pH_Nat2018_10C             = params_Nat2018_NMR_tbl.pH(indices_10C);

Construct = construct_name_Nat2018_all(indices_10C);
pH        = pH_Nat2018_10C;

Nat2018_10C_tbl = cell2table(params_Nat2018_10C_to_save(2:end,:));
Nat2018_10C_tbl.Properties.VariableNames = params_Nat2018_10C_to_save(1,:);
Nat2018_10C_tbl = addvars(Nat2018_10C_tbl,pH,'Before','pGS_pH');
Nat2018_10C_tbl = addvars(Nat2018_10C_tbl,Construct,'Before','pH');

indices_CGC_8p4_10C  = [find(ismember(Construct,{'hpTG-CGC'})&ismember(pH,8.4)),find(ismember(Construct,{'hpUG-CGC'})&ismember(pH,8.4))];
data_CGC_8p4_10C_tbl = Nat2018_10C_tbl(indices_CGC_8p4_10C,:);

%% Figure S5 - Exchange parameters, DNA CGC vs. RNA CGC, pH 8.4, populations

pT_vals_8p4 = [data_CGC_8p4_10C_tbl.pT_pH(1,1);data_CGC_8p4_10C_tbl.pT_pH(2,1)];
pT_errs_8p4 = [data_CGC_8p4_10C_tbl.pT_pH_err(1,1);data_CGC_8p4_10C_tbl.pT_pH_err(2,1)];

pA_vals_8p4 = [data_CGC_8p4_10C_tbl.pA_pH(1,1);data_CGC_8p4_10C_tbl.pA_pH(2,1)];
pA_errs_8p4 = [data_CGC_8p4_10C_tbl.pA_pH_err(1,1);data_CGC_8p4_10C_tbl.pA_pH_err(2,1)];

bar_vals = [pT_vals_8p4.*100,pA_vals_8p4.*100];
bar_errs = [pT_errs_8p4.*100,pA_errs_8p4.*100];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p4;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p4;
hb(2).LineStyle = 'none';

axis([0.5 2.5 0.0 0.5]);

set(gca, 'FontSize', 28, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'Ytick',y_ticks,'TickLength',[0.015,0.015]);

ylabel('population [%]', 'FontSize', 36);

legend({'Tautomer', 'Anion'});

propedit;

%% Figure S5 - Exchange parameters, DNA CGC vs. RNA CGC, pH 8.4, kex

ES1_vals_8p4 = [data_CGC_8p4_10C_tbl.kES1_pH(1,1);data_CGC_8p4_10C_tbl.kES1_pH(2,1)];
ES1_errs_8p4 = [data_CGC_8p4_10C_tbl.kES1_pH_err(1,1);data_CGC_8p4_10C_tbl.kES1_pH_err(2,1)];

ES2_vals_8p4 = [data_CGC_8p4_10C_tbl.kES2_pH(1,1);data_CGC_8p4_10C_tbl.kES2_pH(2,1)];
ES2_errs_8p4 = [data_CGC_8p4_10C_tbl.kES2_pH_err(1,1);data_CGC_8p4_10C_tbl.kES2_pH_err(2,1)];

bar_vals = [ES1_vals_8p4,ES2_vals_8p4];
bar_errs = [ES1_errs_8p4,ES2_errs_8p4];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p4;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p4;
hb(2).LineStyle = 'none';

%axis([0.5 2.5 0.0 0.5]);

set(gca, 'FontSize', 28, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'TickLength',[0.015,0.015]);

ylabel('k_{ex} [s^{-1}]', 'FontSize', 36);

%legend({'Tautomer', 'Anion'});

propedit;

%% Figure S5 - Exchange parameters, DNA CGC vs. RNA CGC, pH 8.4, k_reverse

ES1_kb_vals_8p4 = [data_CGC_8p4_10C_tbl.kES1tGS_pH(1,1);data_CGC_8p4_10C_tbl.kES1tGS_pH(2,1)];
ES1_kb_errs_8p4 = [data_CGC_8p4_10C_tbl.kES1tGS_pH_err(1,1);data_CGC_8p4_10C_tbl.kES1tGS_pH_err(2,1)];

ES2_kb_vals_8p4 = [data_CGC_8p4_10C_tbl.kES2tGS_pH(1,1);data_CGC_8p4_10C_tbl.kES2tGS_pH(2,1)];
ES2_kb_errs_8p4 = [data_CGC_8p4_10C_tbl.kES2tGS_pH_err(1,1);data_CGC_8p4_10C_tbl.kES2tGS_pH_err(2,1)];

bar_vals = [ES1_kb_vals_8p4,ES2_kb_vals_8p4];
bar_errs = [ES1_kb_errs_8p4,ES2_kb_errs_8p4];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p4;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p4;
hb(2).LineStyle = 'none';

%axis([0.5 2.5 0.0 0.5]);

set(gca, 'FontSize', 28, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'TickLength',[0.015,0.015]);

ylabel('k_{reverse} [s^{-1}]', 'FontSize', 36);

%legend({'Tautomer', 'Anion'});

propedit;


%% Figure S5 - Exchange parameters, DNA CGC vs. RNA CGC, pH 8.4, k_forward

ES1_kf_vals_8p4 = [data_CGC_8p4_10C_tbl.kGStES1_pH(1,1);data_CGC_8p4_10C_tbl.kGStES1_pH(2,1)];
ES1_kf_errs_8p4 = [data_CGC_8p4_10C_tbl.kGStES1_pH_err(1,1);data_CGC_8p4_10C_tbl.kGStES1_pH_err(2,1)];

ES2_kf_vals_8p4 = [data_CGC_8p4_10C_tbl.kGStES2_pH(1,1);data_CGC_8p4_10C_tbl.kGStES2_pH(2,1)];
ES2_kf_errs_8p4 = [data_CGC_8p4_10C_tbl.kGStES2_pH_err(1,1);data_CGC_8p4_10C_tbl.kGStES2_pH_err(2,1)];

bar_vals = [ES1_kf_vals_8p4,ES2_kf_vals_8p4];
bar_errs = [ES1_kf_errs_8p4,ES2_kf_errs_8p4];

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

hb = bar(bar_vals);

hold on;

ngroups = size(bar_vals, 1);
nbars = size(bar_vals, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_vals(:,i), bar_errs(:,i), '.','Color',[0 0 0],'LineStyle','none','LineWidth',error_line_width);
end

hold off;

hb(1).FaceColor = 'w';
hb(1).EdgeColor = color_8p4;
hb(1).LineWidth = 2.5;

hb(2).FaceColor = color_8p4;
hb(2).LineStyle = 'none';

%axis([0.5 2.5 0.0 0.5]);

set(gca, 'FontSize', 28, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',2,...
    'TickLength',[0.015,0.015]);

ylabel('k_{forward} [s^{-1}]', 'FontSize', 36);

%legend({'Tautomer', 'Anion'});

propedit;


%% pH interpolation for flux calculations

% Use this script to load raw data from the parent excel file, extract
% kinetic rates and plot figures used for pH extrapolation. Fit ln(k) vs.
% pH to get the slopes and intercepts for each extrapolation.


%% Figure S3 - extrapolation of ln(k) vs. pH

color_ES1 = [0/255,112/255,192/255];
color_ES2 = [0/255,176/255,80/255];

color_ES1tES2 = color_ES1;
color_ES2tES1 = color_ES2;

color_fit_ES1 = [47/255,82/255,143/255];
color_fit_ES2 = [0/255,115/255,80/255];

markersize = 18;

%% Figure S3, left (DNA CGC)

indices_CGC_10C  = find(ismember(Nat2018_10C_tbl.Construct,{'hpTG-CGC'}));

pH_CGC_10C       = Nat2018_10C_tbl.pH(indices_CGC_10C);

kf_ES1_CGC_10C     = Nat2018_10C_tbl.kGStES1_pH(indices_CGC_10C);
kf_ES1_CGC_10C_err = Nat2018_10C_tbl.kGStES1_pH_err(indices_CGC_10C);

kb_ES1_CGC_10C     = Nat2018_10C_tbl.kES1tGS_pH(indices_CGC_10C);
kb_ES1_CGC_10C_err = Nat2018_10C_tbl.kES1tGS_pH_err(indices_CGC_10C);

kf_ES2_CGC_10C     = Nat2018_10C_tbl.kGStES2_pH(indices_CGC_10C);
kf_ES2_CGC_10C_err = Nat2018_10C_tbl.kGStES2_pH_err(indices_CGC_10C);

kb_ES2_CGC_10C     = Nat2018_10C_tbl.kES2tGS_pH(indices_CGC_10C);
kb_ES2_CGC_10C_err = Nat2018_10C_tbl.kES2tGS_pH_err(indices_CGC_10C);

% look at ln(k) for linear fit of ln(k) vs. pH
kf_ES1_CGC_10C_ln     = log(kf_ES1_CGC_10C);
kf_ES1_CGC_10C_ln_err = kf_ES1_CGC_10C_err./kf_ES1_CGC_10C;

kb_ES1_CGC_10C_ln     = log(kb_ES1_CGC_10C);
kb_ES1_CGC_10C_ln_err = kb_ES1_CGC_10C_err./kb_ES1_CGC_10C;

kf_ES2_CGC_10C_ln     = log(kf_ES2_CGC_10C);
kf_ES2_CGC_10C_ln_err = kf_ES2_CGC_10C_err./kf_ES2_CGC_10C;

kb_ES2_CGC_10C_ln     = log(kb_ES2_CGC_10C);
kb_ES2_CGC_10C_ln_err = kb_ES2_CGC_10C_err./kb_ES2_CGC_10C;

% create a matrix of ln(k) values for plotting and fitting
mat_CGC_for_fit_10C = [kf_ES1_CGC_10C_ln, kb_ES1_CGC_10C_ln, kf_ES2_CGC_10C_ln, kb_ES2_CGC_10C_ln];

X = pH_CGC_10C; % x values for fit
init_params = [1 1]; % initial parameters for fitting slope and intercept

x_fit = 7:0.1:10;
y_fit_10C = zeros(length(x_fit),length(mat_CGC_for_fit_10C));

linear_fit_stats_CGC_rates_10C     = zeros(1,length(mat_CGC_for_fit_10C));

slopes_CGC_rates_10C         = zeros(2,length(mat_CGC_for_fit_10C));
intercepts_CGC_rates_10C     = zeros(2,length(mat_CGC_for_fit_10C));

%assume now that kf(ES1) and kb(ES1) are constant because the tautomer
%should not change with pH. Also, it seems that the kf(ES2) changes with
%pH, but the kb(ES2) remains constant

k_cnst = [1, 1, 0, 1]; % logical variable to fit selected ln(k)

for i = 1:length(mat_CGC_for_fit_10C)
    
    if ~k_cnst(i)
        Y = mat_CGC_for_fit_10C(:,i);
        
        f1 = fittype('b+a*x');
        [fit1,gof1,fitinfo1] = fit(X,Y,f1,'StartPoint',init_params); % fit is the function used by cftool
        coeffs1  = coeffvalues(fit1);
        slope1   = coeffs1(1,1);
        inter1   = coeffs1(1,2);
        CI1      = confint(fit1,.95);
        slope1_err = slope1*gof1.sse;
        inter1_err = inter1*gof1.sse;
        
        Rsq1 = gof1.rsquare;
        
        linear_fit_stats_CGC_rates_10C(:,i)     = Rsq1;
        
        y_fit_10C(:,i) = inter1 + slope1.*x_fit';
        
        slopes_CGC_rates_10C(1,i) = slope1;
        slopes_CGC_rates_10C(2,i) = slope1_err;
        intercepts_CGC_rates_10C(1,i) = inter1;
        intercepts_CGC_rates_10C(2,i) = inter1_err;
        
    end
    
end


figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

errorbar(pH_CGC_10C,kf_ES1_CGC_10C_ln,kf_ES1_CGC_10C_ln_err,'o','MarkerEdgeColor',color_ES1,...
    'MarkerFaceColor',color_ES1,'MarkerSize',markersize,'Color',color_ES1,'CapSize',10,'Linewidth',2);
hold on;

errorbar(pH_CGC_10C,kb_ES1_CGC_10C_ln,kb_ES1_CGC_10C_ln_err,'o','MarkerEdgeColor',color_ES1,...
    'MarkerSize',markersize,'Color',color_ES1,'CapSize',10,'Linewidth',2);

errorbar(pH_CGC_10C,kf_ES2_CGC_10C_ln,kf_ES2_CGC_10C_ln_err,'o','MarkerEdgeColor',color_ES2,...
    'MarkerFaceColor',color_ES2,'MarkerSize',markersize,'Color',color_ES2,'CapSize',10,'Linewidth',2);

errorbar(pH_CGC_10C,kb_ES2_CGC_10C_ln,kb_ES2_CGC_10C_ln_err,'o','MarkerEdgeColor',color_ES2,...
    'MarkerSize',markersize,'Color',color_ES2,'CapSize',10,'Linewidth',2);

plot(x_fit,y_fit_10C(:,~k_cnst),'-','Color',color_fit_ES2,'LineWidth',3);

hold off;

axis([7.0 10.0 -2 10]);

set(gca, 'FontSize', 24, 'XminorTick', 'on', 'YminorTick', 'on','LineWidth',2,...
    'TickLength',[0.015,0.015]);

xlabel('pH', 'FontSize', 36);
ylabel('ln(k [s^{-1}])', 'FontSize', 36);

legend({'kf ES1','kb ES1','kf ES2','kb ES2'});

propedit;

%% Figure S3, right (DNA GGC)

indices_GGC_10C  = find(ismember(Nat2018_10C_tbl.Construct,{'hpTG-GGC'}));

pH_GGC_10C       = Nat2018_10C_tbl.pH(indices_GGC_10C);

k_ES2tES1_GGC_10C     = Nat2018_10C_tbl.kES2tES1_pH(indices_GGC_10C);
k_ES2tES1_GGC_10C_err = Nat2018_10C_tbl.kES2tES1_pH_err(indices_GGC_10C);

k_ES1tES2_GGC_10C     = Nat2018_10C_tbl.kES1tES2_pH(indices_GGC_10C);
k_ES1tES2_GGC_10C_err = Nat2018_10C_tbl.kES1tES2_pH_err(indices_GGC_10C);

% look at ln(k) for linear fit of ln(k) vs. pH
k_ES2tES1_GGC_10C_ln = log(k_ES2tES1_GGC_10C);
k_ES2tES1_GGC_10C_ln_err = k_ES2tES1_GGC_10C_err./k_ES2tES1_GGC_10C;

k_ES1tES2_GGC_10C_ln = log(k_ES1tES2_GGC_10C);
k_ES1tES2_GGC_10C_ln_err = k_ES1tES2_GGC_10C_err./k_ES1tES2_GGC_10C;

% create a matrix of ln(k) values for plotting and fitting
mat_GGC_for_fit_10C = [k_ES2tES1_GGC_10C_ln, k_ES1tES2_GGC_10C_ln];

X = pH_GGC_10C;
init_params = [1 1];

x_fit = 7:0.1:10;
y_fit_10C = zeros(length(x_fit),length(mat_GGC_for_fit_10C));

linear_fit_stats_GGC_kminor_10C     = zeros(1,length(mat_GGC_for_fit_10C));

slopes_GGC_kminor_10C         = zeros(2,length(mat_GGC_for_fit_10C));
intercepts_GGC_kminor_10C     = zeros(2,length(mat_GGC_for_fit_10C));

% For kES1-ES2 and KES2-ES1 assume they change with pH

k_cnst = [0, 0]; % logical variable to fit selected ln(k)

for i = 1:length(mat_GGC_for_fit_10C)
    
    if ~k_cnst(i)
        Y = mat_GGC_for_fit_10C(:,i);
        
        f1 = fittype('b+a*x');
        [fit1,gof1,fitinfo1] = fit(X,Y,f1,'StartPoint',init_params); % fit is the function used by cftool
        coeffs1  = coeffvalues(fit1);
        slope1   = coeffs1(1,1);
        inter1   = coeffs1(1,2);
        %CI1      = confint(fit1,.95);
        slope1_err = slope1*gof1.sse;
        inter1_err = inter1*gof1.sse;
        
        Rsq1 = gof1.rsquare;
        
        linear_fit_stats_GGC_kminor_10C(:,i)     = Rsq1;
        
        y_fit_10C(:,i) = inter1 + slope1.*x_fit';
        
        slopes_GGC_kminor_10C(1,i) = slope1;
        slopes_GGC_kminor_10C(2,i) = slope1_err;
        intercepts_GGC_kminor_10C(1,i) = inter1;
        intercepts_GGC_kminor_10C(2,i) = inter1_err;
        
    end
    
end


figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

errorbar(pH_GGC_10C,k_ES2tES1_GGC_10C_ln,k_ES2tES1_GGC_10C_ln_err,'^','MarkerEdgeColor',color_ES2tES1,...
    'MarkerSize',markersize,'Color',color_ES2tES1,'CapSize',10,'Linewidth',2);
hold on;

errorbar(pH_GGC_10C,k_ES1tES2_GGC_10C_ln,k_ES1tES2_GGC_10C_ln_err,'^','MarkerEdgeColor',color_ES1tES2,...
    'MarkerSize',markersize,'Color',color_ES1tES2,'CapSize',10,'Linewidth',2);

plot(x_fit,y_fit_10C(:,1),'-','Color',color_fit_ES2,'LineWidth',3);
plot(x_fit,y_fit_10C(:,2),'-','Color',color_fit_ES1,'LineWidth',3);

hold off;

axis([7.0 10.0 -2 10]);

set(gca, 'FontSize', 24, 'XminorTick', 'on', 'YminorTick', 'on','LineWidth',2,...
    'TickLength',[0.015,0.015]);

xlabel('pH', 'FontSize', 36);
ylabel('ln(k [s^{-1}])', 'FontSize', 36);

legend({'kES2ES1','kES1ES2'});

propedit;


%% fitting 25C data (for extrapolation of hybrid data)

indices_25C = find(ismember(temperatures_Nat2018_all,298.2));

params_Nat2018_25C         = params_Nat2018_NMR_tbl(indices_25C,[25:36,48:59]);
params_Nat2018_25C_to_save = [params_titles;table2cell(params_Nat2018_25C)];
pH_Nat2018_25C             = params_Nat2018_NMR_tbl.pH(indices_25C);

Construct = construct_name_Nat2018_all(indices_25C);
pH        = pH_Nat2018_25C;

Nat2018_25C_tbl = cell2table(params_Nat2018_25C_to_save(2:end,:));
Nat2018_25C_tbl.Properties.VariableNames = params_Nat2018_25C_to_save(1,:);
Nat2018_25C_tbl = addvars(Nat2018_25C_tbl,pH,'Before','pGS_pH');
Nat2018_25C_tbl = addvars(Nat2018_25C_tbl,Construct,'Before','pH');

%% 25C ln(k) vs. pH, DNA CGC

indices_CGC_25C  = find(ismember(Nat2018_25C_tbl.Construct,{'hpTG-CGC'}));

pH_CGC_25C       = Nat2018_25C_tbl.pH(indices_CGC_25C);

kf_ES1_CGC_25C     = Nat2018_25C_tbl.kGStES1_pH(indices_CGC_25C);
kf_ES1_CGC_25C_err = Nat2018_25C_tbl.kGStES1_pH_err(indices_CGC_25C);

kb_ES1_CGC_25C     = Nat2018_25C_tbl.kES1tGS_pH(indices_CGC_25C);
kb_ES1_CGC_25C_err = Nat2018_25C_tbl.kES1tGS_pH_err(indices_CGC_25C);

kf_ES2_CGC_25C     = Nat2018_25C_tbl.kGStES2_pH(indices_CGC_25C);
kf_ES2_CGC_25C_err = Nat2018_25C_tbl.kGStES2_pH_err(indices_CGC_25C);

kb_ES2_CGC_25C     = Nat2018_25C_tbl.kES2tGS_pH(indices_CGC_25C);
kb_ES2_CGC_25C_err = Nat2018_25C_tbl.kES2tGS_pH_err(indices_CGC_25C);

% look at ln(k) for linear fit of ln(k) vs. pH
kf_ES1_CGC_25C_ln     = log(kf_ES1_CGC_25C);
kf_ES1_CGC_25C_ln_err = kf_ES1_CGC_25C_err./kf_ES1_CGC_25C;

kb_ES1_CGC_25C_ln     = log(kb_ES1_CGC_25C);
kb_ES1_CGC_25C_ln_err = kb_ES1_CGC_25C_err./kb_ES1_CGC_25C;

kf_ES2_CGC_25C_ln     = log(kf_ES2_CGC_25C);
kf_ES2_CGC_25C_ln_err = kf_ES2_CGC_25C_err./kf_ES2_CGC_25C;

kb_ES2_CGC_25C_ln     = log(kb_ES2_CGC_25C);
kb_ES2_CGC_25C_ln_err = kb_ES2_CGC_25C_err./kb_ES2_CGC_25C;

% create a matrix of ln(k) values for plotting and fitting
mat_CGC_for_fit_25C = [kf_ES1_CGC_25C_ln, kb_ES1_CGC_25C_ln, kf_ES2_CGC_25C_ln, kb_ES2_CGC_25C_ln];

X = pH_CGC_25C; % x values for fit
init_params = [1 1]; % initial parameters for fitting slope and intercept

x_fit = 7:0.1:10;
y_fit_25C = zeros(length(x_fit),length(mat_CGC_for_fit_25C));

linear_fit_stats_CGC_rates_25C     = zeros(1,length(mat_CGC_for_fit_25C));

slopes_CGC_rates_25C         = zeros(2,length(mat_CGC_for_fit_25C));
intercepts_CGC_rates_25C     = zeros(2,length(mat_CGC_for_fit_25C));

%assume now that kf(ES1) and kb(ES1) are constant because the tautomer
%should not change with pH. Also, it seems that the kf(ES2) changes with
%pH, but the kb(ES2) remains constant

k_cnst = [1, 1, 0, 1]; % logical variable to fit selected ln(k)

for i = 1:length(mat_CGC_for_fit_25C)
    
    if ~k_cnst(i)
        Y = mat_CGC_for_fit_25C(:,i);
        
        f1 = fittype('b+a*x');
        [fit1,gof1,fitinfo1] = fit(X,Y,f1,'StartPoint',init_params); % fit is the function used by cftool
        coeffs1  = coeffvalues(fit1);
        slope1   = coeffs1(1,1);
        inter1   = coeffs1(1,2);
        CI1      = confint(fit1,.95);
        slope1_err = slope1*gof1.sse;
        inter1_err = inter1*gof1.sse;
        
        Rsq1 = gof1.rsquare;
        
        linear_fit_stats_CGC_rates_25C(:,i)     = Rsq1;
        
        y_fit_25C(:,i) = inter1 + slope1.*x_fit';
        
        slopes_CGC_rates_25C(1,i) = slope1;
        slopes_CGC_rates_25C(2,i) = slope1_err;
        intercepts_CGC_rates_25C(1,i) = inter1;
        intercepts_CGC_rates_25C(2,i) = inter1_err;
        
    end
    
end


%% 25C ln(k) vs. pH, DNA GGC

indices_GGC_25C  = find(ismember(Nat2018_25C_tbl.Construct,{'hpTG-GGC'}));

pH_GGC_25C       = Nat2018_25C_tbl.pH(indices_GGC_25C);

k_ES2tES1_GGC_25C     = Nat2018_25C_tbl.kES2tES1_pH(indices_GGC_25C);
k_ES2tES1_GGC_25C_err = Nat2018_25C_tbl.kES2tES1_pH_err(indices_GGC_25C);

k_ES1tES2_GGC_25C     = Nat2018_25C_tbl.kES1tES2_pH(indices_GGC_25C);
k_ES1tES2_GGC_25C_err = Nat2018_25C_tbl.kES1tES2_pH_err(indices_GGC_25C);

% look at ln(k) for linear fit of ln(k) vs. pH
k_ES2tES1_GGC_25C_ln = log(k_ES2tES1_GGC_25C);
k_ES2tES1_GGC_25C_ln_err = k_ES2tES1_GGC_25C_err./k_ES2tES1_GGC_25C;

k_ES1tES2_GGC_25C_ln = log(k_ES1tES2_GGC_25C);
k_ES1tES2_GGC_25C_ln_err = k_ES1tES2_GGC_25C_err./k_ES1tES2_GGC_25C;

% create a matrix of ln(k) values for plotting and fitting
mat_GGC_for_fit_25C = [k_ES2tES1_GGC_25C_ln, k_ES1tES2_GGC_25C_ln];

X = pH_GGC_25C;
init_params = [1 1];

x_fit = 7:0.1:10;
y_fit_25C = zeros(length(x_fit),length(mat_GGC_for_fit_25C));

linear_fit_stats_GGC_kminor_25C     = zeros(1,length(mat_GGC_for_fit_25C));

slopes_GGC_kminor_25C         = zeros(2,length(mat_GGC_for_fit_25C));
intercepts_GGC_kminor_25C     = zeros(2,length(mat_GGC_for_fit_25C));

% For kES1-ES2 and KES2-ES1 assume they change with pH

k_cnst = [0, 0]; % logical variable to fit selected ln(k)

for i = 1:length(mat_GGC_for_fit_25C)
    
    if ~k_cnst(i)
        Y = mat_GGC_for_fit_25C(:,i);
        
        f1 = fittype('b+a*x');
        [fit1,gof1,fitinfo1] = fit(X,Y,f1,'StartPoint',init_params); % fit is the function used by cftool
        coeffs1  = coeffvalues(fit1);
        slope1   = coeffs1(1,1);
        inter1   = coeffs1(1,2);
        %CI1      = confint(fit1,.95);
        slope1_err = slope1*gof1.sse;
        inter1_err = inter1*gof1.sse;
        
        Rsq1 = gof1.rsquare;
        
        linear_fit_stats_GGC_kminor_25C(:,i)     = Rsq1;
        
        y_fit_25C(:,i) = inter1 + slope1.*x_fit';
        
        slopes_GGC_kminor_25C(1,i) = slope1;
        slopes_GGC_kminor_25C(2,i) = slope1_err;
        intercepts_GGC_kminor_25C(1,i) = inter1;
        intercepts_GGC_kminor_25C(2,i) = inter1_err;
        
    end
    
end

%% pH extrapolation for flux calculations, 10C (DNA, RNA)

params_for_flux_calc_10C = data_CGC_8p4_10C_tbl;

pH         = params_for_flux_calc_10C.pH;

pAnion     = 100.*params_for_flux_calc_10C.pA_pH;
pAnion_err = 100.*params_for_flux_calc_10C.pA_pH_err;

pTaut      = 100.*params_for_flux_calc_10C.pT_pH;
pTaut_err  = 100.*params_for_flux_calc_10C.pT_pH_err;

pGS        = 100.*params_for_flux_calc_10C.pGS_pH;
pGS_err    = 100.*params_for_flux_calc_10C.pGS_pH_err;

pKas       = pH - log10(pAnion./(100-pAnion));
pKas_err   = log(10)*sqrt((pAnion_err).^2+(pGS_err).^2);

kGStES1      = params_for_flux_calc_10C.kGStES1_pH;
kGStES1_err  = params_for_flux_calc_10C.kGStES1_pH_err;
kES1tGS      = params_for_flux_calc_10C.kES1tGS_pH;
kES1tGS_err  = params_for_flux_calc_10C.kES1tGS_pH_err;

kGStES2      = params_for_flux_calc_10C.kGStES2_pH;
kGStES2_err  = params_for_flux_calc_10C.kGStES2_pH_err;
kES2tGS      = params_for_flux_calc_10C.kES2tGS_pH;
kES2tGS_err  = params_for_flux_calc_10C.kES2tGS_pH_err;

kES2tES1     = params_for_flux_calc_10C.kES2tES1_pH;
kES2tES1_err = params_for_flux_calc_10C.kES2tES1_pH_err;
kES1tES2     = params_for_flux_calc_10C.kES1tES2_pH;
kES1tES2_err = params_for_flux_calc_10C.kES1tES2_pH_err;

pH_for_calc = 7.4;
params_7p4_for_flux_calc_10C = pH_interpolation(pH_for_calc, pH, pKas, pKas_err,...
    pAnion, pAnion_err, pTaut, pTaut_err, pGS, pGS_err,...
    kGStES1, kGStES1_err, kES1tGS, kES1tGS_err, kGStES2, kGStES2_err, kES2tGS, kES2tGS_err,...
    kES2tES1, kES2tES1_err, kES1tES2, kES1tES2_err,slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

params_7p4_10C_to_save = [params_titles;num2cell(params_7p4_for_flux_calc_10C)];

Construct = params_for_flux_calc_10C.Construct;
pH        = repmat(pH_for_calc,height(params_for_flux_calc_10C),1);

params_7p4_10C_tbl = cell2table(params_7p4_10C_to_save(2:end,:));
params_7p4_10C_tbl.Properties.VariableNames = params_7p4_10C_to_save(1,:);
params_7p4_10C_tbl = addvars(params_7p4_10C_tbl,pH,'Before','pGS_pH');
params_7p4_10C_tbl = addvars(params_7p4_10C_tbl,Construct,'Before','pH');


pH          = params_for_flux_calc_10C.pH;
pH_for_calc = 6.9;
params_6p9_for_flux_calc_10C = pH_interpolation(pH_for_calc, pH, pKas, pKas_err,...
    pAnion, pAnion_err, pTaut, pTaut_err, pGS, pGS_err,...
    kGStES1, kGStES1_err, kES1tGS, kES1tGS_err, kGStES2, kGStES2_err, kES2tGS, kES2tGS_err,...
    kES2tES1, kES2tES1_err, kES1tES2, kES1tES2_err,slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

params_6p9_10C_to_save = [params_titles;num2cell(params_6p9_for_flux_calc_10C)];

Construct = params_for_flux_calc_10C.Construct;
pH        = repmat(pH_for_calc,height(params_for_flux_calc_10C),1);

params_6p9_10C_tbl = cell2table(params_6p9_10C_to_save(2:end,:));
params_6p9_10C_tbl.Properties.VariableNames = params_6p9_10C_to_save(1,:);
params_6p9_10C_tbl = addvars(params_6p9_10C_tbl,pH,'Before','pGS_pH');
params_6p9_10C_tbl = addvars(params_6p9_10C_tbl,Construct,'Before','pH');

%% pH extrapolation for flux calculations, 25C (hybrids)


params_for_flux_calc_25C = Hybrid_8p0_25C_tbl;

pH         = params_for_flux_calc_25C.pH;

pAnion     = 100.*params_for_flux_calc_25C.pA_pH;
pAnion_err = 100.*params_for_flux_calc_25C.pA_pH_err;

pTaut      = 100.*params_for_flux_calc_25C.pT_pH;
pTaut_err  = 100.*params_for_flux_calc_25C.pT_pH_err;

pGS        = 100.*params_for_flux_calc_25C.pGS_pH;
pGS_err    = 100.*params_for_flux_calc_25C.pGS_pH_err;

pKas       = pH - log10(pAnion./(100-pAnion));
pKas_err   = log(10)*sqrt((pAnion_err).^2+(pGS_err).^2);

kGStES1      = params_for_flux_calc_25C.kGStES1_pH;
kGStES1_err  = params_for_flux_calc_25C.kGStES1_pH_err;
kES1tGS      = params_for_flux_calc_25C.kES1tGS_pH;
kES1tGS_err  = params_for_flux_calc_25C.kES1tGS_pH_err;

kGStES2      = params_for_flux_calc_25C.kGStES2_pH;
kGStES2_err  = params_for_flux_calc_25C.kGStES2_pH_err;
kES2tGS      = params_for_flux_calc_25C.kES2tGS_pH;
kES2tGS_err  = params_for_flux_calc_25C.kES2tGS_pH_err;

kES2tES1     = params_for_flux_calc_25C.kES2tES1_pH;
kES2tES1_err = params_for_flux_calc_25C.kES2tES1_pH_err;
kES1tES2     = params_for_flux_calc_25C.kES1tES2_pH;
kES1tES2_err = params_for_flux_calc_25C.kES1tES2_pH_err;

pH_for_calc = 7.4;
params_7p4_for_flux_calc_25C = pH_interpolation(pH_for_calc, pH, pKas, pKas_err,...
    pAnion, pAnion_err, pTaut, pTaut_err, pGS, pGS_err,...
    kGStES1, kGStES1_err, kES1tGS, kES1tGS_err, kGStES2, kGStES2_err, kES2tGS, kES2tGS_err,...
    kES2tES1, kES2tES1_err, kES1tES2, kES1tES2_err,slopes_CGC_rates_25C, slopes_GGC_kminor_25C);

params_7p4_25C_to_save = [params_titles;num2cell(params_7p4_for_flux_calc_25C)];

Construct = params_for_flux_calc_25C.Construct;
pH        = repmat(pH_for_calc,height(params_for_flux_calc_25C),1);

params_7p4_25C_tbl = cell2table(params_7p4_25C_to_save(2:end,:));
params_7p4_25C_tbl.Properties.VariableNames = params_7p4_25C_to_save(1,:);
params_7p4_25C_tbl = addvars(params_7p4_25C_tbl,pH,'Before','pGS_pH');
params_7p4_25C_tbl = addvars(params_7p4_25C_tbl,Construct,'Before','pH');


pH          = params_for_flux_calc_25C.pH;
pH_for_calc = 6.9;
params_6p9_for_flux_calc_25C = pH_interpolation(pH_for_calc, pH, pKas, pKas_err,...
    pAnion, pAnion_err, pTaut, pTaut_err, pGS, pGS_err,...
    kGStES1, kGStES1_err, kES1tGS, kES1tGS_err, kGStES2, kGStES2_err, kES2tGS, kES2tGS_err,...
    kES2tES1, kES2tES1_err, kES1tES2, kES1tES2_err,slopes_CGC_rates_25C, slopes_GGC_kminor_25C);

params_6p9_25C_to_save = [params_titles;num2cell(params_6p9_for_flux_calc_25C)];

Construct = params_for_flux_calc_25C.Construct;
pH        = repmat(pH_for_calc,height(params_for_flux_calc_25C),1);

params_6p9_25C_tbl = cell2table(params_6p9_25C_to_save(2:end,:));
params_6p9_25C_tbl.Properties.VariableNames = params_6p9_25C_to_save(1,:);
params_6p9_25C_tbl = addvars(params_6p9_25C_tbl,pH,'Before','pGS_pH');
params_6p9_25C_tbl = addvars(params_6p9_25C_tbl,Construct,'Before','pH');

%% Create parameter tables for flux calcualtions

params_for_flux_calc_all  = [data_CGC_8p4_10C_tbl;Hybrid_8p0_25C_tbl;...
    params_7p4_10C_tbl;Hybrid_25C_tbl(1,:);params_7p4_25C_tbl(2,:);...
    params_6p9_10C_tbl;params_6p9_25C_tbl];
params_for_flux_calc_CGC_DNA  = data_CGC_8p4_10C_tbl(1,:);
params_for_flux_calc_CGC_RNA  = data_CGC_8p4_10C_tbl(2,:);

% save the tables as .csv files for python to read for simulations
subfolder       ='params_for_python_flux';
subfolder_path  = sprintf('%s/%s',parent_path,subfolder);

writetable(params_for_flux_calc_all,sprintf('%s/params_for_flux_calc_all.csv',subfolder_path));
% open this file using python to simulate popuplations at equilibrim using
% the MIS scheme (including k2). Results are saved as csv files under
% Flux-Simulations-Pol-Epsilon/k2s_all/

writetable(params_for_flux_calc_CGC_DNA,sprintf('%s/params_for_flux_calc_CGC_DNA.csv',subfolder_path));
writetable(params_for_flux_calc_CGC_RNA,sprintf('%s/params_for_flux_calc_CGC_RNA.csv',subfolder_path));
% open these files using python to simulate popuplations at equilibrim using
% the MIS scheme (including k1 and k2). Results are saved as csv files under
% Flux-Simulations-Pol-Epsilon/k2s_k_1s_all/

%% Figure 4 and Figure S7 - flux vs. k2

path_flux_sims   = 'Flux-Simulations-Pol-Epsilon/k2s_all/';
parent_path_sims = sprintf('%s/%s',parent_path,path_flux_sims);

filenames = {'CGC_DNA_8p4','CGC_RNA_8p4','dGrU_pH_8p0','dTrG_pH_8p0',...
    'CGC_DNA_7p4','CGC_RNA_7p4','dGrU_pH_7p4','dTrG_pH_7p4',...
    'CGC_DNA_6p9','CGC_RNA_6p9','dGrU_pH_6p9','dTrG_pH_6p9'}';

select_filenames_pH  = filenames(1:4,1);
select_filenames_7p4 = filenames(5:8,1);
select_filenames_6p9 = filenames(9:12,1);

colors_lines = {[0/255 0/255 255/255];[255/255 0/255 0/255];...
    [255/255 0/255 255/255];[170/255 0/255 200/255]};

pol_epsilon_k2 = 268; %s-1
pol_beta_k2    = 1365; %s-1
pol_T7_k2      = 660; %s-1

% create k2 vector
length_k2 = 115;
k2_first = 0.001;
k2_multip = 1.235;

k2s = zeros(1,length_k2);
k2s(1,1) = k2_first;

for i = 2:length_k2
    k2s(1,i) = k2_multip*k2s(1,i-1);
end


%% Figure 4 - plot fA vs k2s - 7p4

fAs_pH = zeros(length(select_filenames_7p4),length(k2s));

for i = 1:length(select_filenames_7p4)
    
    path_data = sprintf('%s%s.csv',parent_path_sims,select_filenames_7p4{i,1});
    rawData = readtable(path_data);
    
    fAs_pH(i,:) = rawData.fA;

end

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

h2 = plot(k2s,(fAs_pH.*100),'LineWidth',3);

set(h2, {'color'}, colors_lines);

hold on;

xline(pol_epsilon_k2,'k-','LineWidth',3);
xline(pol_beta_k2,'k-','LineWidth',3);
xline(pol_T7_k2,'k-','LineWidth',3);

hold off;

set(gca, 'FontSize', 18, 'LineWidth', 2, 'TickLength',[0.015,0.015],...
      'YminorTick', 'on', 'XScale', 'log');
  
axis([1E-3 2E7 0.0 90]);

legend(select_filenames_7p4); 

propedit;


%% Figure S7A - plot fA vs k2s - 6p9

fAs_6p9 = zeros(length(select_filenames_6p9),length(k2s));



for i = 1:length(select_filenames_6p9)
    
    path_data = sprintf('%s%s.csv',parent_path_sims,select_filenames_6p9{i,1});
    rawData = readtable(path_data);
    
    fAs_6p9(i,:) = rawData.fA;

end

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

h3 = plot(k2s,(fAs_6p9.*100),'LineWidth',3);

set(h3, {'color'}, colors_lines);

hold on;

xline(pol_epsilon_k2,'k-','LineWidth',3);
xline(pol_beta_k2,'k-','LineWidth',3);
xline(pol_T7_k2,'k-','LineWidth',3);

hold off;

set(gca, 'FontSize', 18, 'LineWidth', 2, 'TickLength',[0.015,0.015],...
      'YminorTick', 'on', 'XScale', 'log');

axis([1E-3 2E7 0.0 90]);

legend(select_filenames_6p9); 

propedit;

%% Figure S7B - plot fA vs k2s - high pHs

fAs_pH = zeros(length(select_filenames_pH),length(k2s));

for i = 1:length(select_filenames_pH)
    
    path_data = sprintf('%s%s.csv',parent_path_sims,select_filenames_pH{i,1});
    rawData = readtable(path_data);
    
    fAs_pH(i,:) = rawData.fA;

end

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

h1 = plot(k2s,(fAs_pH.*100),'LineWidth',3);

set(h1, {'color'}, colors_lines);

hold on;

xline(pol_epsilon_k2,'k-','LineWidth',3);
xline(pol_beta_k2,'k-','LineWidth',3);
xline(pol_T7_k2,'k-','LineWidth',3);

hold off;

set(gca, 'FontSize', 18, 'LineWidth', 2, 'TickLength',[0.015,0.015],...
      'YminorTick', 'on', 'XScale', 'log');
  
axis([1E-3 2E7 0.0 90]);

legend(select_filenames_pH); 

propedit;

%% Figure S6 - flux vs. k2 and k_1

path_flux_sims   = 'Flux-Simulations-Pol-Epsilon/k2s_k_1s_all/';
parent_path_sims = sprintf('%s/%s',parent_path,path_flux_sims);

path_CGC_DNA  = 'sims_CGC_DNA';
path_CGC_RNA  = 'sims_CGC_RNA';

pol_epsilon_k_1_dTTP = 70000; %s-1

% create k_1 vector
length_k_1 = 130; %130
k_1_first = 25.0; %1.0
k_1_multip = 1.125; %1.15

k_1s = zeros(1,length_k_1);
k_1s(1,1) = k_1_first;

for i = 2:length_k_1
    k_1s(1,i) = k_1_multip*k_1s(1,i-1);
end


%% Figure S6B left - heatmap CGC DNA

path_data = sprintf('%s%s',parent_path_sims,path_CGC_DNA);
formatSpec = '%.1f';

fA_CGC_DNA = zeros(length_k2,length_k_1);

for i=1:length_k_1
    
    k_1_str = num2str(k_1s(1,i),formatSpec);
    filename = sprintf('%s_%s.csv',path_CGC_DNA,k_1_str);
    filepath = sprintf('%s/%s',path_data,filename);
    
    rawData = readtable(filepath);
    
    fA_CGC_DNA(:,i) = rawData.fA;
    
end

z_min = 0;
z_max = 100;

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

h1 = imagesc(k2s,k_1s,(fA_CGC_DNA.*100)');
hold on;
plot(pol_epsilon_k2,pol_epsilon_k_1_dTTP,'ko','MarkerFaceColor','k');
hold off;

caxis([z_min, z_max]);

set(gca, 'FontSize', 18, 'LineWidth', 2, 'TickLength',[0.015,0.015],...
      'YDir', 'normal', 'YScale', 'log', 'XScale', 'log');
%set(gca, 'YDir', 'normal', 'YScale', 'log', 'XScale', 'log');

colormap(jet(256));
colorbar;

propedit;

%% Figure S6B right - heatmap CGC RNA

path_data = sprintf('%s%s',parent_path_sims,path_CGC_RNA);
formatSpec = '%.1f';

fA_CGC_RNA = zeros(length_k2,length_k_1);

for i=1:length_k_1
    
    k_1_str = num2str(k_1s(1,i),formatSpec);
    filename = sprintf('%s_%s.csv',path_CGC_RNA,k_1_str);
    filepath = sprintf('%s/%s',path_data,filename);
    
    rawData = readtable(filepath);
    
    fA_CGC_RNA(:,i) = rawData.fA;
    
end

z_min = 0;
z_max = 100;

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

h2 = imagesc(k2s,k_1s,(fA_CGC_RNA.*100)');
hold on;
plot(pol_epsilon_k2,pol_epsilon_k_1_dTTP,'ko','MarkerFaceColor','k');
hold off;

caxis([z_min, z_max]);

set(gca, 'FontSize', 18, 'LineWidth', 2, 'TickLength',[0.015,0.015],...
      'YDir', 'normal', 'YScale', 'log', 'XScale', 'log');
%set(gca, 'YDir', 'normal', 'YScale', 'log', 'XScale', 'log');

colormap(jet(256));
colorbar;

propedit;

%% Figure S1 - NMR Spectra

parent_path_NMR = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-NMR-Spectra';
path_1D_2D_imino = '2023-1D-2D-imino-RawData';

path_7p4 = sprintf('%s/%s/pH_7p4',parent_path_NMR,path_1D_2D_imino);
path_7p8 = sprintf('%s/%s/pH_7p8',parent_path_NMR,path_1D_2D_imino);

color_dTrG_pH_7p4 = [200/255,40/255,100/255];
color_dGrU_pH_7p4 = [200/255,40/255,100/255];
color_dTrG_pH_7p8 = [46/255,49/255,146/255];
color_dGrU_pH_7p8 = [46/255,49/255,146/255];

linewidth_1D   = 3;
linewidth_plot = 2;

%% Figure S1A

filename = 'dTrG_pH7p4_1D.mat';
load(sprintf('%s/%s',path_7p4,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

plot(dTrG_pH7p4_1D{1,1},dTrG_pH7p4_1D{1,2},'-',...
    'Color',color_dTrG_pH_7p4,'LineWidth',linewidth_1D);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', 'YTick', [],...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 -0.5E5 6E5]);

propedit;

filename = 'dGrU_pH7p4_1D.mat';
load(sprintf('%s/%s',path_7p4,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

plot(dGrU_pH7p4_1D{1,1},dGrU_pH7p4_1D{1,2},'-',...
    'Color',color_dGrU_pH_7p4,'LineWidth',linewidth_1D);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', 'YTick', [],...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 -1E5 11E5]);

propedit;


%% Figure S1B

filename = 'dTrG_pH7p4_2D.mat';
load(sprintf('%s/%s',path_7p4,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

contour(dTrG_pH7p4_2D{1,1},dTrG_pH7p4_2D{1,2},...
    dTrG_pH7p4_2D{1,3}',dTrG_pH7p4_2D{1,4},...
    'Color',color_dTrG_pH_7p4);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', ...
    'Ydir', 'rev', 'YTick', 140:5:165, 'YminorTick', 'off',...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 138 162]);

propedit;


filename = 'dGrU_pH7p4_2D.mat';
load(sprintf('%s/%s',path_7p4,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

contour(dGrU_pH7p4_2D{1,1},dGrU_pH7p4_2D{1,2},...
    dGrU_pH7p4_2D{1,3}',dGrU_pH7p4_2D{1,4},...
    'Color',color_dGrU_pH_7p4);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', ...
    'Ydir', 'rev', 'YTick', 140:5:165, 'YminorTick', 'off',...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 138 162]);

propedit;


%% Figure S1C

filename = 'dTrG_pH7p8_1D.mat';
load(sprintf('%s/%s',path_7p8,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

plot(dTrG_pH7p8_1D{1,1},dTrG_pH7p8_1D{1,2},'-',...
    'Color',color_dTrG_pH_7p8,'LineWidth',linewidth_1D);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', 'YTick', [],...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 -2E5 5E6]);

propedit;

filename = 'dGrU_pH7p8_1D.mat';
load(sprintf('%s/%s',path_7p8,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

plot(dGrU_pH7p8_1D{1,1},dGrU_pH7p8_1D{1,2},'-',...
    'Color',color_dGrU_pH_7p8,'LineWidth',linewidth_1D);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', 'YTick', [],...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 -5E5 9.5E6]);

propedit;

%% Figure S1D

filename = 'dTrG_pH7p8_2D.mat';
load(sprintf('%s/%s',path_7p8,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

contour(dTrG_pH7p8_2D{1,1},dTrG_pH7p8_2D{1,2},...
    dTrG_pH7p8_2D{1,3}',dTrG_pH7p8_2D{1,4},...
    'Color',color_dTrG_pH_7p8);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', ...
    'Ydir', 'rev', 'YTick', 140:5:165, 'YminorTick', 'off',...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 138 162]);

propedit;


filename = 'dGrU_pH7p8_2D.mat';
load(sprintf('%s/%s',path_7p8,filename));

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

contour(dGrU_pH7p8_2D{1,1},dGrU_pH7p8_2D{1,2},...
    dGrU_pH7p8_2D{1,3}',dGrU_pH7p8_2D{1,4},...
    'Color',color_dGrU_pH_7p8);

set(gca, 'FontSize', 28, 'Xdir', 'rev', 'XminorTick', 'off', ...
    'Ydir', 'rev', 'YTick', 140:5:165, 'YminorTick', 'off',...
    'LineWidth', linewidth_plot, 'TickLength', [0.015,0.015]);

axis([9 14.5 138 162]);

propedit;

%% Figrue 3 - NTP NMR titrations and fitting titrations to get NTP pKa

color_rUTP = [46/255,49/255,146/255];
color_dTTP = [180/255,50/255,120/255];

parent_path_titrations = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2023/2023-RNA-DNA-Hybrid/2023-Manuscript-Figures/2023-For-Data-Deposition/2023-NMR-Spectra/2023-NMR-titrations';
file_name_pHs          = 'dTTP-rUTP-actual-pH.xlsx';

filename_dTTP = 'cs_dTTP.mat';
load(sprintf('%s/%s',parent_path_titrations,filename_dTTP));

pH_NMRtube_dTTP = readtable(sprintf('%s/%s',parent_path_titrations,file_name_pHs),'Sheet','dTTP');

pH_dTTP    = pH_NMRtube_dTTP.actualPH;
d_obs_dTTP = cs(:,4);

sigmoid_init_params = [0 10 12];% L {func(1)} & H {func(2)} & pKa {func(3)}

sigmoid_fit_func    = @(func,x_data) func(1) + (func(2)-func(1))./(1+10.^(x_data-func(3)))+ 1000*(func(1) < 0)^2 + 1000*(func(3) < 0)^2; % func = L + ((H-L)/(1+10^(x-pka)))

fault = 0;

x_data = pH_dTTP;
y_data = d_obs_dTTP;

try
    [sigmoid_fitting_params_dTTP(1,:), r_sigmoid, J_sigmoid, covB_sigmoid] = nlinfit(x_data,y_data,sigmoid_fit_func,sigmoid_init_params); %non-linear fit
    sigmoid_Rsq_dTTP(1,:)       = 1 - sum(r_sigmoid.^2) / sum((y_data - mean(y_data)).^2);
    sigmoid_covMats_dTTP{1,1}   = sqrt(covB_sigmoid);
    sigmoid_error_pka_dTTP(1,:) = sqrt(covB_sigmoid(3,3));
    sigmoid_error_H1_dTTP(1,:)   = sqrt(covB_sigmoid(2,2));
    sigmoid_error_L1_dTTP(1,:)   = sqrt(covB_sigmoid(1,1));
catch %#ok<CTCH>
    fault=fault+1;
    sigmoid_fitting_params_dTTP(1,:)=[0,0,0];
    sigmoid_Rsq_dTTP(1,:)=[0];
end



filename_rUTP = 'cs_rUTP.mat';
load(sprintf('%s/%s',parent_path_titrations,filename_rUTP));

pH_NMRtube_rUTP = readtable(sprintf('%s/%s',parent_path_titrations,file_name_pHs),'Sheet','rUTP');

pH_rUTP    = pH_NMRtube_rUTP.ActualPH;
d_obs_rUTP = cs(:,4);

fault = 0;

x_data = pH_rUTP;
y_data = d_obs_rUTP;

try
    [sigmoid_fitting_params_rUTP(1,:), r_sigmoid, J_sigmoid, covB_sigmoid] = nlinfit(x_data,y_data,sigmoid_fit_func,sigmoid_init_params); %non-linear fit
    sigmoid_Rsq_rUTP(1,:)       = 1 - sum(r_sigmoid.^2) / sum((y_data - mean(y_data)).^2);
    sigmoid_covMats_rUTP{1,1}   = sqrt(covB_sigmoid);
    sigmoid_error_pka_rUTP(1,:) = sqrt(covB_sigmoid(3,3));
    sigmoid_error_H1_rUTP(1,:)   = sqrt(covB_sigmoid(2,2));
    sigmoid_error_L1_rUTP(1,:)   = sqrt(covB_sigmoid(1,1));
catch %#ok<CTCH>
    fault=fault+1;
    sigmoid_fitting_params_rUTP(1,:)=[0,0,0];
    sigmoid_Rsq_rUTP(1,:)=[0];
end

pKa_dTTP = [sigmoid_fitting_params_dTTP(1,3),sigmoid_error_pka_dTTP];
pKa_rUTP = [sigmoid_fitting_params_rUTP(1,3),sigmoid_error_pka_rUTP];

x_fits = 0:0.1:14;
y_fits_dTTP = sigmoid_fit_func(sigmoid_fitting_params_dTTP(1,:),x_fits);
y_fits_rUTP = sigmoid_fit_func(sigmoid_fitting_params_rUTP(1,:),x_fits);


figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

plot(pH_rUTP,d_obs_rUTP,'o','MarkerEdgeColor',color_rUTP,'MarkerFaceColor',color_rUTP,'MarkerSize',14);

hold on;

plot(pH_dTTP,d_obs_dTTP,'o','MarkerEdgeColor',color_dTTP,'MarkerFaceColor',color_dTTP,'MarkerSize',14);

plot(x_fits,y_fits_rUTP,'-','Color',color_rUTP,'LineWidth',2.5);
plot(x_fits,y_fits_dTTP,'-','Color',color_dTTP,'LineWidth',2.5);

hold off;

set(gca, 'FontSize', 28, 'XminorTick', 'on', 'YminorTick', 'on','LineWidth',2,...
    'TickLength',[0.015,0.015]);

axis([6.4 14 -0.5 10]);

xlabel('pH', 'FontSize', 36);
ylabel('Chemical Shift difference [ppm]', 'FontSize', 36);

propedit;