%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MYOPIC LABOR UNION     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE MANAGEMENT AND HOUSEKEEPING      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close;
timer_tot=tic;

% Which parts to run and print
run_optss=1;
print_optss=1;
run_init=1;
print_init=1;
run_findpath=1;
print_findpath=1;
run_optpath=1;
print_optpath=1;
print_startend=1;
run_compiler=0;

% Which type of graphs to print
print_graph_type = '-depsc';
%print_graph_type = '-dpdf';

% Create folder structure
iteration=4; % Where to save files
path = sprintf('optimal_tax_test%d',iteration);
mkdir(path);
path_docs=path;
path_graphs=[path,'/Graphs'];
path_graphs_rel = 'Graphs';
mkdir(path_docs);
mkdir(path_graphs);
copyfile('opttax_noc_parameters.m',[path_docs,'/','parameters.m']);
diary(sprintf('diary%d',iteration));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEADY STATE TAX RATE UNION       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------');
disp('STEADY STATE TAX RATE');
disp('--------------------------------');

run opttax_noc_parameters.m;
run opttax_noc_steadystate.m;

tau_ss_LU = taul;
g_a_ss_LU = sum(g.*dz,2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEADY STATE NO TAX                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run opttax_noc_parameters.m;
taul = 0;
gg_sparse=sparse(ones(I*J,1)./(I*J*dz*da_stacked)); %Uniform guess
run opttax_noc_steadystate_no_LU.m;

g_a_ss_NT = sum(g.*dz,2);
g_z = g'*da_tilde;
disp('Done!');

% Save figure comparing SS wealth distributions
H=figure;
set(gcf,'Visible','on');
plot(a,g_a_ss_NT,'-b',a,g_a_ss_LU,'-r');
legend('No tax',sprintf('Union tax = %d',tau_ss_LU));
title('Steady state wealth distributions');
xlabel('Wealth a');
ylabel('Density');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_wealthdistss%d',iteration)]);
delete(H);

% Save productivity distribution
H=figure;
set(gcf,'Visible','off');
plot(z,g_z,'-b');
title('Steady state productivity distribution');
xlabel('Productivity z');
ylabel('Density');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_productivitydistss%d',iteration)]);
delete(H);




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL DISTRIBUTION                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------');
disp('FINDING INITIAL DISTRIBUTION');
disp('--------------------------------');

run opttax_noc_parameters.m;
taul = 0;
gg_sparse=sparse(ones(I*J,1)./(I*J*dz*da_stacked)); %Uniform guess
run opttax_noc_steadystate_no_LU.m;

disp('Done!');

% Option 4: Start from steady state distribution with no tax, scaled down by "factor"
at=factor*a;
ind=zeros(I,1);
rest=zeros(I,1);
for i=1:I
    if at(i)>=min(a)
        ind(i)=max(find(a-at(i)<=0));
        rest(i)=(at(i)-a(ind(i)))/(a(ind(i)+1)-a(ind(i)));
    else
        ind(i)=1;
        rest(i)=0;
    end
    
    if ind(i)>=I
        ind(i)=da-1;
        rest(i)=1;
        printf('Warning: index>=da')
    end
        
    if rest(i)>1
        rest(i)=1;
        printf('Warning: rest>1')
    end
end
g0=zeros(I,J);
for i=1:I
    g0(ind(i),:)=g0(ind(i),:)+(1-rest(i))*g(i,:);
    g0(ind(i)+1,:)=g0(ind(i)+1,:)+rest(i)*g(i,:);
end

% MARGINALS AND SAVE
% Reshape and save distribution
gg0 = g0(:);
gg0 = gg0/(gg0'*da_stacked*dz);
g_a01 = sum(g.*dz,2);
g_a02 = sum(g0.*dz,2);
gg_initial = gg0;
gg_guess = gg;

% SAVE FIGURES
% Create figures
H=figure;
set(gcf,'Visible','off');
plot(a,g_a02,'-b');
xlim([0 amax]);
title('Initial wealth distribution');
xlabel('Wealth a');
ylabel('Density');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_initialwealthdistss%d',iteration)]);
delete(H);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSITION NO TAX           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('----------------------------');
disp('PLOTTING TRANSITION RESULTS');

% RUN AND SAVE RESULTS FOR LAISSEZ-FAIRE (NO TAX)  
disp('NO TAX');
% Find transition path
taul_bar = 0; 
taul0 = 0;
gal = 0;
taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
taulvec_lf = taulvec;
disp('Iteration    Max wage update    Max abs. ES');
[welfare_mat,gg,w_t,pi_t,~,ls_it] = opttax_noc_transition_no_LU(gg_guess,gg_initial,taulvec);
% Save results
run opttax_measures.m;
profits_total_vec_lf = profits_total_vec;
profit_share_vec_lf = profit_share_vec;
frac_active_vec_lf = frac_active_vec;
frac_constrained_vec_lf = frac_constrained_vec;
frac_output_constrained_vec_lf = frac_output_constrained_vec;
wealth_tot_vec_lf = wealth_tot_vec;
capital_vec_lf = capital_vec;
labor_vec_lf = labor_vec;
output_vec_lf = output_vec;
tfp_vec_lf = tfp_vec;
utility_vec_lf = utility_vec;
frac_entrepreneurs_lf = frac_entrepreneurs;
frac_output_of_unconstrained_vec_lf = frac_output_of_unconstrained_vec;
wage_vec_lf = w_t;
welfare_workers_lf=welfare_mat(:,2);
welfare_entrepreneurs_lf=welfare_mat(:,3);
welfare_lf = (1-paretoweight_workers)*welfare_mat(:,3)+paretoweight_workers*welfare_mat(:,2);
DY = dy;
pi_lf = pi_t;
save tax_iteration_results_lf.mat profits_total_vec_lf profit_share_vec_lf frac_active_vec_lf frac_constrained_vec_lf frac_output_constrained_vec_lf wealth_tot_vec_lf ...
    capital_vec_lf labor_vec_lf output_vec_lf tfp_vec_lf utility_vec_lf frac_entrepreneurs_lf frac_output_of_unconstrained_vec_lf wage_vec_lf ...
    welfare_workers_lf welfare_entrepreneurs_lf welfare_lf taulvec_lf DY pi_lf; 
disp('Done!');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSITION UNION TAX           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RUN AND SAVE RESULTS FOR OPTIMAL TAX SCHEDULE
disp('UNION TAX SCHEDULE');
% Find transition path
disp('Iteration    Max wage update    Max abs. ES   Max abs. pi diff');
[welfare_mat,gg,w_t,pi_t,taulvec,flag] = opttax_noc_transition(gg_initial);
% Save results
run opttax_measures.m;
profits_total_vec_opt = profits_total_vec;
profit_share_vec_opt = profit_share_vec;
frac_active_vec_opt = frac_active_vec;
frac_constrained_vec_opt = frac_constrained_vec;
frac_output_constrained_vec_opt = frac_output_constrained_vec;
wealth_tot_vec_opt = wealth_tot_vec;
capital_vec_opt = capital_vec;
labor_vec_opt = labor_vec;
output_vec_opt = output_vec;
tfp_vec_opt = tfp_vec;
utility_vec_opt = utility_vec;
frac_entrepreneurs_opt = frac_entrepreneurs;
frac_output_of_unconstrained_vec_opt = frac_output_of_unconstrained_vec;
wage_vec_opt = w_t;
welfare_workers_opt=welfare_mat(:,2);
welfare_entrepreneurs_opt=welfare_mat(:,3);
welfare_opt = (1-paretoweight_workers)*welfare_mat(:,3)+paretoweight_workers*welfare_mat(:,2);
pi_opt = pi_t;
taulvec_opt = taulvec;
save tax_iteration_results_opt.mat profits_total_vec_opt profit_share_vec_opt frac_active_vec_opt frac_constrained_vec_opt frac_output_constrained_vec_opt wealth_tot_vec_opt ...
    capital_vec_opt labor_vec_opt output_vec_opt tfp_vec_opt utility_vec_opt frac_entrepreneurs_opt frac_output_of_unconstrained_vec_opt wage_vec_opt...
    welfare_workers_opt welfare_entrepreneurs_opt welfare_opt taulvec_opt pi_opt; 
disp('Done!');



% CREATE PLOTS COMPARING RESULTS 
load tax_iteration_results_lf.mat;
load tax_iteration_results_opt.mat;

% Figures with optimal tax as deviation from no tax
H=figure;
set(gcf,'Visible','off');
set(gca,'FontSize',10);
%
subplot(3,2,1); plot(timevec,(labor_vec_opt-labor_vec_lf)./abs(labor_vec_lf),'-g',timevec,zeros(N,1),'--b');
xlim([0 T]); xlabel('Years');
title('(a) Labor supply, l');
%
subplot(3,2,2); plot(timevec,(wealth_tot_vec_opt-wealth_tot_vec_lf)./abs(wealth_tot_vec_lf),'-g',timevec,zeros(N,1),'--b');
xlim([0 T]); xlabel('Years');
title('(b) Entrepreneurial wealth, x');
%
subplot(3,2,3); plot(timevec,(output_vec_opt./labor_vec_opt-output_vec_lf./labor_vec_lf)./(output_vec_lf./labor_vec_lf),'-r',...
    timevec,(wage_vec_opt-wage_vec_lf)./abs(wage_vec_lf),'-g',...
    timevec,zeros(N,1),'--b');
xlim([0 T]); xlabel('Years');
title('(c) Wage, w, and labor productivity, y/l');
%
subplot(3,2,4); plot(timevec,(tfp_vec_opt-tfp_vec_lf)./abs(tfp_vec_lf),'-g',timevec,zeros(N,1),'--b');
xlim([0 T]); xlabel('Years');
title('(d) Total factor productivity, Z');
%
subplot(3,2,5); plot(timevec,(output_vec_opt-output_vec_lf)./abs(output_vec_lf),'-g',timevec,zeros(N,1),'--b');
xlim([0 T]); xlabel('Years');
title('(e) Income, y');
%
subplot(3,2,6); plot(timevec,(utility_vec_opt-utility_vec_lf)./abs(utility_vec_lf),'-g',timevec,zeros(N,1),'--b');
xlim([0 T]); xlabel('Years');
title('(f) Worker period utility, u(c,l)');
%
print(print_graph_type,[path_graphs,'/',sprintf('fig_deviationoptlf%d',iteration)]);
delete(H);


% Figures with both optimal tax and no tax
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,taulvec_opt,'-b');
xlabel('Time');
title(sprintf('Union tax schedule'));
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_taxschedule%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,pi_opt,'-b',timevec,pi_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Labor employed by unconstrained entrepreneurs');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_pischedule%d',iteration)]);
delete(H);


H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,wage_vec_opt,'-b',timevec,wage_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Wage');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_wage%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,profits_total_vec_opt,'-b',timevec,profits_total_vec_lf,'-r');
legend('Optimal tax','No tax');
xlabel('Time');
title('Total profits');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_profitstot%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,profit_share_vec_opt,'-b',timevec,profit_share_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Profits/(profits + labor income)');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_profitshare%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_entrepreneurs_opt,'-b',timevec,frac_entrepreneurs_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Fraction entrepreneurs');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracentrepreneurs%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_constrained_vec_opt,'-b',timevec,frac_constrained_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Fraction entrepreneurs constrained');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracconstrained%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_output_constrained_vec_opt,'-b',timevec,frac_output_constrained_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Fraction output from constrained entrepreneurs');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracoutputconstrained%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,labor_vec_opt,'-b',timevec,labor_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Aggregate labor supply');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_labor%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,capital_vec_opt,'-b',timevec,capital_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Aggregate capital use');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_capital%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,output_vec_opt,'-b',timevec,output_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('GDP');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_gdp%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,tfp_vec_opt,'-b',timevec,tfp_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('TFP');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_tfp%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_output_of_unconstrained_vec_opt,'-b',timevec,frac_output_of_unconstrained_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('GDP/Potential GDP');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracgdpofpotential%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,wealth_tot_vec_opt,'-b',timevec,wealth_tot_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Aggregate wealth');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_aggwealth%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,utility_vec_opt,'-b',timevec,utility_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Period utility');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_periodutility%d',iteration)]);
delete(H);
   
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,(1-taulvec_opt).*wage_vec_opt.*labor_vec_opt,'-b',timevec,(1-taulvec_lf).*wage_vec_lf.*labor_vec_lf,'-r');
legend('Union tax','No tax');
xlabel('Time');
title('Worker period consumption');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_periodconsumption%d',iteration)]);
delete(H);

save plotting_vectors.mat timevec labor_vec_opt labor_vec_lf wealth_tot_vec_opt wealth_tot_vec_lf output_vec_opt output_vec_lf wage_vec_opt wage_vec_lf ...
tfp_vec_opt tfp_vec_lf utility_vec_opt utility_vec_lf taulvec_opt profits_total_vec_opt profits_total_vec_lf profit_share_vec_opt profit_share_vec_lf ...
frac_entrepreneurs_opt frac_entrepreneurs_lf frac_constrained_vec_opt frac_constrained_vec_lf frac_output_constrained_vec_opt frac_output_constrained_vec_lf ...
capital_vec_opt capital_vec_lf frac_output_of_unconstrained_vec_opt frac_output_of_unconstrained_vec_lf pi_opt pi_lf DY;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE LATEX FILE        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVE PLOTS TO LATEX FILE
if print_optpath==1
    
load tax_iteration_results_lf.mat;
load tax_iteration_results_opt.mat;


copyfile('output.tex','texfile.tex');
fid=fopen('texfile.tex','a');

% Table comparing welfare of various tax schedules

%load tax_iteration_results_optflat.mat;

% fprintf(fid,'\\begin{table}[p] \n');
% fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | c | c | c } \n \\hline \n');
% fprintf(fid,' & No tax & Optimal policy & Optimal SS tax & Best flat policy \\\\ \n');
% fprintf(fid,'\\hline \n');
% fprintf(fid,'$\\tau_0$ & %7.5f & %7.5f & %7.5f & %7.5f  \\\\ \n',taulvec_lf(1),taulvec_opt(1),taulvec_optss(1),taulvec_optflat(1));
% fprintf(fid,'$\\bar{\\tau}$ & %7.5f & %7.5f & %7.5f & %7.5f \\\\ \n',taulvec_lf(end),taulvec_opt(end),taulvec_optss(end),taulvec_optflat(end));
% fprintf(fid,'Half life & - & %7.5f & - & -  \\\\ \n',halflife_opt);
% fprintf(fid,'Welfare & %7.5f & %7.5f & %7.5f & %7.5f \\\\ \n',welfare_lf,welfare_opt,welfare_optss,welfare_optflat);
% fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n');
% fprintf(fid,'\\end{table} \n \n');

fprintf(fid,'\\section{Tables with  Results} \n '); 

%fprintf(fid,'\\begin{table}[p]\n');
fprintf(fid,'\\begin{table}[h!]\n');
fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | c } \n \\hline \n');
fprintf(fid,' & No tax & Union tax \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'$\\tau_0$ & %7.5f & %7.5f \\\\ \n',taulvec_lf(1),taulvec_opt(1));
fprintf(fid,'$\\bar{\\tau}$ & %7.5f & %7.5f \\\\ \n',taulvec_lf(end),taulvec_opt(end));
fprintf(fid,'Welfare (weighted) & %7.5f & %7.5f \\\\ \n',welfare_lf,welfare_opt);
fprintf(fid,'Welfare workers & %7.5f & %7.5f \\\\ \n',welfare_workers_lf,welfare_workers_opt);
fprintf(fid,'Welfare entrepreneurs & %7.5f & %7.5f \\\\ \n',welfare_entrepreneurs_lf,welfare_entrepreneurs_opt);
fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n');
fprintf(fid,'\\end{table} \n \n');


%fprintf(fid,'Optimal flat tax rate = %5.3f \n',opttaul_flat);

rho_ave = (pop_share/discrate_workers+(1-pop_share)/rho)^(-1);

fprintf(fid,'\\begin{table}[h!]\n');
fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | c | c } \n \\hline \n');
fprintf(fid,' Experiment & Total welfare & Worker welfare & Entrepreneur welfare \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'Optimal policy & %7.5f & %7.5f & %7.5f \\\\ \n',exp(rho_ave*(welfare_opt-welfare_lf))-1,exp(discrate_workers*(welfare_workers_opt-welfare_workers_lf))-1,exp(rho*(welfare_entrepreneurs_opt-welfare_entrepreneurs_lf))-1);
fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n');
fprintf(fid,'\\end{table} \n \n');

fprintf(fid,'\\clearpage \n \n');

if print_startend==1
    
% PRINT PARAMETERS ETC.
run opttax_noc_parameters.m;
% Print some choices
fprintf(fid,'\\section{Parameters and functional forms} \n '); 
fprintf(fid,'\\subsection{Functional forms etc.} \n '); 
fprintf(fid,'\\begin{itemize} \n');
fprintf(fid,'\\item Occupational choice: No \n \\item Workers save: No \n \\item Decreasing returns to scale: Yes \n');
fprintf(fid,'\\item Productivity process: Ornstein-Uhlenbeck, $ d\\log{(z)}= -\\nu\\log{(z)} dt + \\sigma dW $ \n');
if preftype==1
    fprintf(fid,'\\item Period utility function: \n $$ u(c,l)=(1-\\gamma)^{-1}(c-\\nu(l))^{1-\\gamma}, \\quad \\nu(l)=(1+1/\\chi)^{-1}l^{1+1/\\chi} $$ \n');
elseif preftype==2
    fprintf(fid,'\\item Period utility function: \n $$ u(c,l)=(1-\\gamma)^{-1}c^{1-\\gamma}-\\nu(l), \\quad \\nu(l)=(1+1/\\chi)^{-1}l^{1+1/\\chi} $$ \n');
end
fprintf(fid,'\\item Production function: $ y = F(z,k,n) = zA((k-f_k)^{+})^{\\alpha}((n-f_n)^{+})^{\\beta} $ \n');
fprintf(fid,'\\item Tax schedule: $ \\tau_l(t) = \\bar{\\tau}_l + e^{-\\gamma t} (\\tau_{l,0}-\\bar{\\tau}_l) $ \n');
fprintf(fid,'\\end{itemize} \n');
% Print parameter values
fprintf(fid,'\\subsection{Parameter values} \n '); 
fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | r } \n \\hline \n');
fprintf(fid,'Pareto weight workers &  & %5.3f \\\\ \n',paretoweight_workers);
fprintf(fid,'Population share of workers & $popshare$ & %5.3f \\\\ \n',pop_share);
fprintf(fid,'Total population & $popmass$ & %5.3f \\\\ \n',pop_mass);
fprintf(fid,'Discount rate entrepreneurs & $\\rho_e$ & %5.3f \\\\ \n',rho);
fprintf(fid,'Discount rate workers & $\\rho_w$ & %5.3f \\\\ \n',discrate_workers);
fprintf(fid,'Relative risk aversion & $\\gamma$ & %5.3f \\\\ \n',gam);
fprintf(fid,'Labor disutility parameter & $\\chi $ & %5.3f \\\\ \n',chi0);
fprintf(fid,'Depreciation rate & $\\delta$ & %5.3f \\\\ \n',d);
fprintf(fid,'Death rate & $\\theta$ & %5.3f \\\\ \n',the);
fprintf(fid,'Fixed cost capital & $f_k$ & %5.3f \\\\ \n',fk);
fprintf(fid,'Fixed cost labor & $f_n$ & %5.3f \\\\ \n',fn);
fprintf(fid,'Financial constraint parameter & $\\lambda$ & %5.3f \\\\ \n',la);
fprintf(fid,'Common TFP parameter & $A$ & %5.3f \\\\ \n',Aprod);
fprintf(fid,'Capital share & $\\alpha$ & %5.3f \\\\ \n',al);
fprintf(fid,'Labor share & $\\beta$ & %5.3f \\\\ \n',be);
fprintf(fid,'Returns to scale & $\\alpha+\\beta$ & %5.3f \\\\ \n',al+be);
fprintf(fid,'Interest rate & $r^*$ & %5.3f \\\\ \n',rstar);
fprintf(fid,'Effect of productivity on effective labor supply & $\\eta$ & %5.3f \\\\ \n',eta);
fprintf(fid,'Productivity drift parameter & $\\nu$ & %5.3f \\\\ \n',nu);
fprintf(fid,'Productivity yearly autocorrelation & $e^{-\\nu}$ & %5.3f \\\\ \n',exp(-nu));
fprintf(fid,'Productivity standard deviation parameter & $\\sigma$ & %5.3f \\\\ \n',sqrt(sig2));
fprintf(fid,'Productivity mean & $\\bar{z}$ & %5.3f \\\\ \n',exp(logzmean + Var/2));
fprintf(fid,'Poisson arrival rate & &  %5.3f \\\\ \n',poissonarrivalrate*poissonshocks);
if poissonshocks==1
    fprintf(fid,'Parameter of Pareto distribution of Poisson shocks & &  %5.3f \\\\ \n',paretoparam);
end
fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n \n');
% Print iteration parameter
fprintf(fid,'\\subsection{Iteration parameters} \n ');
fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | r } \n \\hline \n');
fprintf(fid,'Number of grid points assets & $I$ & %5.3f \\\\ \n',I);
fprintf(fid,'Number of grid points productivity & $J$ & %5.3f \\\\ \n',J);
fprintf(fid,'Number of grid points time & $N$ & %5.3f \\\\ \n',N);
fprintf(fid,'Number of time periods & $T$ & %5.3f \\\\ \n',T);
fprintf(fid,'Max assets & $a_{max}$ & %5.3f \\\\ \n',amax);
fprintf(fid,'Mean wealth relative to steady state &  & %5.3f \\\\ \n',contract);
fprintf(fid,'Contraction of initial distribution & $factor$ & %5.3f \\\\ \n',factor);
fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n \n');

fprintf(fid,'\\clearpage \n \n');

% Section: figures
fprintf(fid,'\\section{Figures} \n');

end

% Print Optimal Tax Rate Figures
if print_optss==1

fprintf(fid,'Union steady state tax rate = %5.3f \n',tau_ss_LU);
% Productivity distribution
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_productivitydistss%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n');
% Wealth distributions
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_wealthdistss%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n');
end

% Print initial wealth distribution
if print_init==1
% Save initial wealth distribution
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_initialwealthdistss%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 
end

fprintf(fid,'\\clearpage \n \n');

fprintf(fid,'\\clearpage \n \n');

% Figures with optimal tax as deviation from no tax
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_deviationoptlf%d',iteration)]);
fprintf(fid,'\\caption{Proportional deviations of optimal tax equilibrium from the laissez-faire equilibrium}\n');
fprintf(fid,'\\end{figure}\n\n'); 

% Figures with both no tax and optimal tax
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_taxschedule%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_pischedule%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_wage%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n');

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_profitstot%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_profitshare%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_fracentrepreneurs%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_fracconstrained%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_fracoutputconstrained%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_labor%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_capital%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_gdp%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_tfp%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_fracgdpofpotential%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_aggwealth%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_periodutility%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 

fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_periodconsumption%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n');

end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT END OF FILE                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if print_startend==1
    
% End latex document
fprintf(fid,'\n \\end{document}');
fclose(fid);

% Copy and move files into appropriate folder
movefile(sprintf('diary%d',iteration),path_docs), %Move diary
movefile('texfile.tex',path_docs); %Move file to correct subfolder
copyfile('main_Labor_Union.m',path_docs);

end

diary off; %Close diary


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPILATION AND FILE MANAGEMENT       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run Latex compiler and move files
if run_compiler==1
    
% Run Latex compiler
setenv('PATH', [getenv('PATH') ':/usr/local/texlive/2016/bin/x86_64-darwin']); %Put latex compiler in path (depends on local system)
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
% setenv('PATH', [getenv('PATH') ':/usr/bin']); %For Nobel server
!pdflatex texfile.tex > Texjunk.junk
% !latex output.tex > Texjunk.junk
% !dvips output.dvi > DVIjunk.junk
% !ps2pdf output.ps > PSjunk.junk

% Delete, move and rename files
delete('texfile.aux');
delete('texfile.log');
delete('texfile.out');
delete('Texjunk.junk');

end



