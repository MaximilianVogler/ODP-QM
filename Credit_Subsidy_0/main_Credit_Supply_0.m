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
% OPTIMAL STEADY STATE CREDIT SUBSIDY RATE         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_optss==1
    
disp('--------------------------------');
disp('OPTIMAL STEADY STATE CREDIT SUBSIDY RATE');
disp('--------------------------------');

% Find optimal credit subsidy by iteration over grid
run opttax_noc_parameters.m;
g_optss_struct = cell(length(tax_vec_ss));
welfare_vec_ss = zeros(length(tax_vec_ss),1);
disp(' Credit Subsidy      Welfare total    Welfare workers    Welfare entrepreneurs');
for ittax=1:length(tax_vec_ss)
    taul = 0;
    credk=tax_vec_ss(ittax);
    gg_sparse=sparse(ones(I*J,1)./(I*J*dz*da_stacked)); %Uniform guess
    run opttax_noc_steadystate.m;
    welfare_tot_ss = (1-paretoweight_workers)*welfare_entrepreneurs+paretoweight_workers*welfare_workers;
    welfare_vec_ss(ittax) = welfare_tot_ss;
    g_optss_struct{ittax} = g;
    fprintf('%5.3f        %7.5f          %7.5f             %7.5f \n',credk,welfare_tot_ss,welfare_workers,welfare_entrepreneurs);
end
disp('Done!');
[maxwelfare_ss,maxwelfare_ss_idx]=max(welfare_vec_ss);
opttaul_ss=tax_vec_ss(maxwelfare_ss_idx);
g_optss = g_optss_struct{maxwelfare_ss_idx};
fprintf('Optimal credit subsidy = %5.3f \n',opttaul_ss);
g_a_opttax = sum(g.*dz,2);


% Find no credit subsidy results
disp('RESULTS WITH NO CREDIT SUBSIDY');
taul = 0;
credk = 0;
gg_sparse=sparse(ones(I*J,1)./(I*J*dz*da_stacked));
run opttax_noc_steadystate.m;
g_a_notax = sum(g.*dz,2);
g_z = g'*da_tilde;
disp('Done!');

% Save figure comparing wealth distributions
H=figure;
set(gcf,'Visible','on');
plot(a,g_a_notax,'-b',a,g_a_opttax,'-r');
legend('No credit subsidy',sprintf('Optimal credit subsidy = %d',opttaul_ss));
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

% Save figure of welfare as function of credit subsidy
H=figure;
set(gcf,'Visible','off');
plot(tax_vec_ss,welfare_vec_ss,'--b');
title('Worker welfare as function of credit subsidy');
xlabel('Credit Subsidy $\varsigma_k$','interpreter','latex');
ylabel('Worker welfare');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_welfaress%d',iteration)]);
delete(H);

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL DISTRIBUTION                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_init==1

disp('--------------------------------');
disp('FINDING INITIAL DISTRIBUTION');

% DISTRIBUTION TO START FROM
run opttax_noc_parameters.m;
taul = 0;
credk = 0;
gg_sparse=sparse(ones(I*J,1)./(I*J*dz*da_stacked)); %Uniform guess
run opttax_noc_steadystate.m;

disp('Done!');

% CONVERT TO INITIAL DISTRIBUTION
% Option 1: Point mass at fraction of mean wealth of original distribution
% agg_wealth = (as.*g(:))'*da_stacked*dz;
% mean_productivity = (zs.*g(:))'*da_stacked*dz;
% g_z = da_tilde'*g;
% [aux,ind1]=min(abs(a-contract*agg_wealth));
% [aux,ind2]=min(abs(z-mean_productivity));
% g0 = zeros(I,J);
% %g0(ind1,ind2)=1/(da_tilde(ind1)*dz); %Point mass at mean productivity
% g0(ind1,:) = g_z/da_tilde(ind1); %Marginal of productivity is unchanged

% Option 2: Contract all wealth point by a given amount from original distribution
% [aux,idx]=min(abs(repmat(a*contract,1,I)-repmat(a',I,1)),[],2);
% g0 = zeros(I,J);
% for jj=1:J
% for ii=1:I
%     g0(idx(ii),jj)=g0(idx(ii),jj)+g(ii,jj);
% end
% end

% Option 3: Start from steady state distribution under the optimal steady
% state credit subsidy
% g0=g_optss;

% Option 4: Start from steady state distribution with no credit subsidy, scaled down by "factor"
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

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMAL TRANSITION FLAT CREDIT SUBSIDY           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matrices for iteration
[taubar_mat_tr,tau0_mat_tr,gal_mat_tr] = meshgrid(taubar_vec_tr,tau0_vec_tr,gal_vec_tr);
[gp_tau0_tr,gp_taubar_tr,gp_gal_tr]=size(taubar_mat_tr);
gp_tot_tr=gp_taubar_tr*gp_tau0_tr*gp_gal_tr;
%halflives_mat = repmat(reshape(halflives_vec,1,1,gp_gal),gp_tau0,gp_taubar,1);

%Saving results
welfare_workers_mat1_tr = zeros(gp_tau0_tr,gp_taubar_tr,gp_gal_tr);
welfare_workers_mat2_tr = zeros(gp_tau0_tr,gp_taubar_tr,gp_gal_tr);
welfare_entrepreneurs_tr = zeros(gp_tau0_tr,gp_taubar_tr,gp_gal_tr);
flag_mat_tr=zeros(gp_tau0_tr,gp_taubar_tr,gp_gal_tr);
maxES_mat_tr=zeros(gp_tau0_tr,gp_taubar_tr,gp_gal_tr);
%dcompletediter_tr=zeros(gp_tau0_tr,gp_taubar_tr,gp_gal_tr);

% Loop over credit subsidy parameters
timer_findpath = tic;
% fprintf('Progress:\n');
% fprintf(['\n' repmat('.',1,gp_tot) '\n\n']);
disp('Iteration    Max wage update    Max ES');
parfor ittax=1:gp_tot_tr
    
    % PARAMETERS
    taul = 0;
    taul_bar_tr = taubar_mat_tr(ittax); 
    taul0_tr = tau0_mat_tr(ittax);
    gal_tr = gal_mat_tr(ittax);
    taulvec_tr = taul_bar_tr + exp(-gal_tr.*timevec).*(taul0_tr-taul_bar_tr);

    % TRANSITION PATH
    [welfare_mat_tr,~,~,~,flag] = opttax_noc_transition(gg_guess,gg_initial,taulvec_tr);
    
    % Save and print progress
    welfare_workers_mat1_tr(ittax)=welfare_mat_tr(:,1);
    welfare_workers_mat2_tr(ittax)=welfare_mat_tr(:,2);
    welfare_entrepreneurs_tr(ittax)=welfare_mat_tr(:,3);
    flag_mat_tr(ittax)=flag; %=1 if no convergence in steady state, =2 if no convergence in transition
    fprintf('Iteration %d/%d done! \n',ittax,gp_tot_tr);
end
disp('Done!');
fprintf('Total time to find optimal credit subsidy schedule = %3.0f minutes and %2.0f seconds \n',floor(toc(timer_findpath)/60),rem(toc(timer_findpath),60))

% Welfare measures
welfaremat_workers_tr = welfare_workers_mat2_tr;
welfaremat_tr = (1-paretoweight_workers)*welfare_entrepreneurs_tr+paretoweight_workers*welfaremat_workers_tr; %Approximation to worker welfare used
welfaremat_relative_tr = ((welfaremat_tr-max(max(max(welfaremat_tr))))/abs(max(max(max(welfaremat_tr)))))*100;
welfaremat_workers_relative_tr = ((welfaremat_workers_tr-max(max(max(welfaremat_workers_tr))))/abs(max(max(max(welfaremat_workers_tr)))))*100;
welfaremat_entrepreneurs_relative_tr = ((welfare_entrepreneurs_tr-max(max(max(welfare_entrepreneurs_tr))))/abs(max(max(max(welfare_entrepreneurs_tr)))))*100;
% Save optimal credit subsidy schedule
idx=find(welfaremat_tr==max(max(max(welfaremat_tr))));
tau0_opt_tr=tau0_mat_tr(idx);
taubar_opt_tr=taubar_mat_tr(idx);
%halflife_opt=halflives_mat(idx);
gal_opt_tr = gal_mat_tr(idx);
% Save optimal flat credit subsidy schedule
[aux,idx]=max(diag(welfaremat_tr(:,:,1)));
opttaul_flat=tau0_vec_tr(idx);

% Save measures
save([path_docs,'/','findpath_results_tr.mat'],'welfaremat_tr','tau0_mat_tr','taubar_mat_tr','tau0_opt_tr','taubar_opt_tr','gal_opt_tr','opttaul_flat');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMAL TRANSITION CREDIT SUBSIDY SCHEDULE       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_findpath==1
    
disp('-------------------------------------------');
disp('FINDING THE OPTIMAL TRANSITION CREDIT SUBSIDY SCHEDULE');

% General parameters
run opttax_noc_parameters.m;

%Matrices for iteration
[taubar_mat,tau0_mat,gal_mat] = meshgrid(taubar_vec,tau0_vec,gal_vec);
[gp_tau0,gp_taubar,gp_gal]=size(taubar_mat);
gp_tot=gp_taubar*gp_tau0*gp_gal;
halflives_mat = repmat(reshape(halflives_vec,1,1,gp_gal),gp_tau0,gp_taubar,1);

%Saving results
welfare_workers_mat1 = zeros(gp_tau0,gp_taubar,gp_gal);
welfare_workers_mat2 = zeros(gp_tau0,gp_taubar,gp_gal);
welfare_entrepreneurs = zeros(gp_tau0,gp_taubar,gp_gal);
flag_mat=zeros(gp_tau0,gp_taubar,gp_gal);
maxES_mat=zeros(gp_tau0,gp_taubar,gp_gal);
dcompletediter=zeros(gp_tau0,gp_taubar,gp_gal);

% Loop over credit subsidy parameters
timer_findpath = tic;
% fprintf('Progress:\n');
% fprintf(['\n' repmat('.',1,gp_tot) '\n\n']);
disp('Iteration    Max wage update    Max ES');
parfor ittax=1:gp_tot

    
    % PARAMETERS
    taul = 0;
    taul_bar = taubar_mat(ittax); 
    taul0 = tau0_mat(ittax);
    gal = gal_mat(ittax);
    taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);

    % TRANSITION PATH
    [welfare_mat,~,~,~,flag] = opttax_noc_transition(gg_guess,gg_initial,taulvec);
    
    % Save and print progress
    welfare_workers_mat1(ittax)=welfare_mat(:,1);
    welfare_workers_mat2(ittax)=welfare_mat(:,2);
    welfare_entrepreneurs(ittax)=welfare_mat(:,3);
    flag_mat(ittax)=flag; %=1 if no convergence in steady state, =2 if no convergence in transition
    %maxES_mat(ittax)=maxESdiff;
    dcompleted(ittax)=1;
    %fprintf('Iteration %d/%d done!',sum(sum(sum(dcompleted))),gp_tot);
    fprintf('Iteration %d/%d done! \n',ittax,gp_tot);
end
disp('Done!');
fprintf('Total time to find optimal credit subsidy schedule = %3.0f minutes and %2.0f seconds \n',floor(toc(timer_findpath)/60),rem(toc(timer_findpath),60))

% Welfare measures
welfaremat_workers = welfare_workers_mat2;
welfaremat = (1-paretoweight_workers)*welfare_entrepreneurs+paretoweight_workers*welfaremat_workers; %Approximation to worker welfare used
welfaremat_relative = ((welfaremat-max(max(max(welfaremat))))/abs(max(max(max(welfaremat)))))*100;
welfaremat_workers_relative = ((welfaremat_workers-max(max(max(welfaremat_workers))))/abs(max(max(max(welfaremat_workers)))))*100;
welfaremat_entrepreneurs_relative = ((welfare_entrepreneurs-max(max(max(welfare_entrepreneurs))))/abs(max(max(max(welfare_entrepreneurs)))))*100;
% Save optimal credit subsidy schedule
idx=find(welfaremat==max(max(max(welfaremat))));
tau0_opt=tau0_mat(idx);
taubar_opt=taubar_mat(idx);
halflife_opt=halflives_mat(idx);
gal_opt = gal_mat(idx);
% Save optimal flat credit subsidy schedule
% if (length(tau0_vec)==length(taubar_vec)) && (tau0_vec==taubar_vec)
%     [aux,idx]=max(diag(welfaremat(:,:,1)));
%     opttaul_flat=tau0_vec(idx);
% end
% Save measures
save([path_docs,'/','findpath_results.mat'],'welfaremat','tau0_mat','taubar_mat','halflives_mat','tau0_opt','taubar_opt','halflife_opt','gal_opt');

% PLOT WELFARE
for ii=1:gp_gal
    H=figure;
    set(gcf,'Visible','off');
    xlim([taubar_vec(1) taubar_vec(end)]);
    ylim([tau0_vec(1) tau0_vec(end)]);
    %Relative welfare: total
    surfc(taubar_mat(:,:,ii),tau0_mat(:,:,ii),welfaremat_relative(:,:,ii));
    xlabel('Long-run credit subsidy');
    ylabel('Initial credit subsidy');
    zlabel('Worker welfare (per cent deviation from optimum)');
    title(sprintf('Welfare by credit subsidy schedule, half life = %d',halflives_vec(ii)));
    set(gca,'FontSize',14);
    print(print_graph_type,[path_graphs,'/',sprintf('fig_welfarepathtax_hl%d_%d',halflives_vec(ii),iteration)]);
    savefig([path_graphs,'/',sprintf('fig_welfarepathtax_hl%d_%d',halflives_vec(ii),iteration)]);
    delete(H);
    %Relative welfare: workers
    H=figure;
    surfc(taubar_mat(:,:,ii),tau0_mat(:,:,ii),welfaremat_workers_relative(:,:,ii));
    xlabel('Long-run credit subsidy');
    ylabel('Initial credit subsidy');
    zlabel('Worker welfare (per cent deviation from optimum)');
    title(sprintf('Welfare by credit subsidy schedule, half life = %d',halflives_vec(ii)));
    set(gca,'FontSize',14);
    print(print_graph_type,[path_graphs,'/',sprintf('fig_welfarepathtax_workers_hl%d_%d',halflives_vec(ii),iteration)]);
    savefig([path_graphs,'/',sprintf('fig_welfarepathtax_workers_hl%d_%d',halflives_vec(ii),iteration)]);
    delete(H);
    %Relative welfare: entrepreneurs
    H=figure;
    surfc(taubar_mat(:,:,ii),tau0_mat(:,:,ii),welfaremat_entrepreneurs_relative(:,:,ii));
    xlabel('Long-run credit subsidy');
    ylabel('Initial credit subsidy');
    zlabel('Entrepreneur welfare (per cent deviation from optimum)');
    title(sprintf('Welfare by credit subsidy schedule, half life = %d',halflives_vec(ii)));
    set(gca,'FontSize',14);
    print(print_graph_type,[path_graphs,'/',sprintf('fig_welfarepathtax_entrepreneurs_hl%d_%d',halflives_vec(ii),iteration)]);
    savefig([path_graphs,'/',sprintf('fig_welfarepathtax_entrepreneurs_hl%d_%d',halflives_vec(ii),iteration)]);
    delete(H);
    %Total welfare: total
    surfc(taubar_mat(:,:,ii),tau0_mat(:,:,ii),welfaremat(:,:,ii));
    xlabel('Long-run credit subsidy');
    ylabel('Initial credit subsidy');
    zlabel('Total welfare');
    title(sprintf('Welfare by credit subsidy schedule, half life = %d',halflives_vec(ii)));
    set(gca,'FontSize',14);
    print(print_graph_type,[path_graphs,'/',sprintf('fig_welfarepathtax_hl_total%d_%d',halflives_vec(ii),iteration)]);
    savefig([path_graphs,'/',sprintf('fig_welfarepathtax_hl_total%d_%d',halflives_vec(ii),iteration)]);
    delete(H);
end

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS FOR OPTIMAL CREDIT SUBSIDY AND NO CREDIT SUBSIDY      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_optpath==1

load([path_docs,'/','findpath_results.mat']);
    
disp('----------------------------');
disp('PLOTTING TRANSITION RESULTS');

% RUN AND SAVE RESULTS FOR LAISSEZ-FAIRE (NO CREDIT SUBSIDY)  
disp('NO CREDIT SUBSIDY');
% Find transition path
taul = 0;
taul_bar = 0; 
taul0 = 0;
gal = 0;
taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
taulvec_lf = taulvec;
disp('Iteration    Max wage update    Max abs. ES');
[welfare_mat,gg,w_t,consumption_workers,~,ls_it] = opttax_noc_transition(gg_guess,gg_initial,taulvec);
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
save tax_iteration_results_lf.mat profits_total_vec_lf profit_share_vec_lf frac_active_vec_lf frac_constrained_vec_lf frac_output_constrained_vec_lf wealth_tot_vec_lf ...
    capital_vec_lf labor_vec_lf output_vec_lf tfp_vec_lf utility_vec_lf frac_entrepreneurs_lf frac_output_of_unconstrained_vec_lf wage_vec_lf ...
    welfare_workers_lf welfare_entrepreneurs_lf welfare_lf taulvec_lf DY; 
disp('Done!');

% RUN AND SAVE RESULTS FOR OPTIMAL CREDIT SUBSIDY SCHEDULE
disp('OPTIMAL CREDIT SUBSIDY SCHEDULE');
% Find transition path
taul = 0;
taul_bar = taubar_opt; 
taul0 = tau0_opt;
gal = gal_opt;
taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
taulvec_opt = taulvec;
disp('Iteration    Max wage update    Max abs. ES');
[welfare_mat,gg,w_t,consumption_workers,~,ls_it] = opttax_noc_transition(gg_guess,gg_initial,taulvec);
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
save tax_iteration_results_opt.mat profits_total_vec_opt profit_share_vec_opt frac_active_vec_opt frac_constrained_vec_opt frac_output_constrained_vec_opt wealth_tot_vec_opt ...
    capital_vec_opt labor_vec_opt output_vec_opt tfp_vec_opt utility_vec_opt frac_entrepreneurs_opt frac_output_of_unconstrained_vec_opt wage_vec_opt...
    welfare_workers_opt welfare_entrepreneurs_opt welfare_opt taulvec_opt halflife_opt; 
disp('Done!');


% RUN AND SAVE RESULTS FOR FIXED CREDIT SUBSIDY AT STEADY STATE OPTIMUM
disp('FIXED CREDIT SUBSIDY: STEADY STATE OPTIMUM');
% Find transition path
taul = 0;
taul_bar = opttaul_ss; 
taul0 = opttaul_ss;
gal = 0;
taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
taulvec_optss = taulvec;
disp('Iteration    Max wage update    Max abs. ES');
[welfare_mat,~,~,~,~] = opttax_noc_transition(gg_guess,gg_initial,taulvec);
% Save results
welfare_workers_optss=welfare_mat(:,2);
welfare_entrepreneurs_optss=welfare_mat(:,3);
welfare_optss = (1-paretoweight_workers)*welfare_mat(:,3)+paretoweight_workers*welfare_mat(:,2);
save tax_iteration_results_optss.mat welfare_workers_optss welfare_entrepreneurs_optss welfare_optss taulvec_optss;
disp('Done!');

% RUN AND SAVE RESULTS FOR OPTIMAL FLAT CREDIT SUBSIDY
disp('FIXED CREDIT SUBSIDY: STEADY STATE OPTIMUM');
% Find transition path
taul = 0;
taul_bar = opttaul_flat; 
taul0 = opttaul_flat;
gal = 0;
taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
taulvec_flat = taulvec;
disp('Iteration    Max wage update    Max abs. ES');
[welfare_mat,~,~,~,~] = opttax_noc_transition(gg_guess,gg_initial,taulvec);
% Save results
welfare_workers_flat=welfare_mat(:,2);
welfare_entrepreneurs_flat=welfare_mat(:,3);
welfare_flat = (1-paretoweight_workers)*welfare_mat(:,3)+paretoweight_workers*welfare_mat(:,2);
save tax_iteration_results_flat.mat welfare_workers_flat welfare_entrepreneurs_flat welfare_flat taulvec_flat;
disp('Done!');


% RUN AND SAVE RESULTS FOR OPTIMAL FLAT CREDIT SUBSIDY
% disp('FIXED CREDIT SUBSIDY: OPTIMAL FLAT CREDIT SUBSIDY');
% % Find transition path
% taul_bar = opttaul_flat; 
% taul0 = opttaul_flat;
% gal = 0;
% taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
% taulvec_optflat = taulvec;
% [welfare_mat,~,~,~,~] = opttax_noc_transition(gg_guess,gg_initial,taulvec);
% % Save results
% welfare_optflat = welfare_mat(:,2);
% save tax_iteration_results_optflat.mat welfare_optflat taulvec_optflat;
% disp('Converged');

% RUN AND SAVE RESULTS FOR constant tau0_opt
taul = 0;
taul_bar = tau0_opt; 
taul0 = tau0_opt;
gal = 0;
taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
taulvec_t0 = taulvec;
disp('Iteration    Max wage update    Max abs. ES');
[welfare_mat,~,~,~,~] = opttax_noc_transition(gg_guess,gg_initial,taulvec_t0);
% Save results
welfare_workers_t0=welfare_mat(:,2);
welfare_entrepreneurs_t0=welfare_mat(:,3);
welfare_t0 = (1-paretoweight_workers)*welfare_mat(:,3)+paretoweight_workers*welfare_mat(:,2);
save tax_iteration_results_t0.mat welfare_workers_t0 welfare_entrepreneurs_t0 welfare_t0 taulvec_t0;
disp('Done!');

% RUN AND SAVE RESULTS FOR constant taubar_opt
taul = 0;
taul_bar = taulvec_opt(end); 
taul0 = taulvec_opt(end);
gal = 0;
taulvec = taul_bar + exp(-gal.*timevec).*(taul0-taul_bar);
taulvec_bar = taulvec;
disp('Iteration    Max wage update    Max abs. ES');
[welfare_mat,~,~,~,~] = opttax_noc_transition(gg_guess,gg_initial,taulvec_bar);
% Save results
welfare_workers_bar=welfare_mat(:,2);
welfare_entrepreneurs_bar=welfare_mat(:,3);
welfare_bar = (1-paretoweight_workers)*welfare_mat(:,3)+paretoweight_workers*welfare_mat(:,2);
save tax_iteration_results_bar.mat welfare_workers_bar welfare_entrepreneurs_bar welfare_bar taulvec_bar;
disp('Done!');

% COMPUTE SS D/Y FOR OPTIMAL TAU_BAR
% run opttax_noc_parameters.m;
% taul=taulvec_opt(end);
% gg_sparse=sparse(ones(I*J,1)./(I*J*dz*da_stacked)); %Uniform guess
% run opttax_noc_steadystate.m;
% GDP = (ddprod(:).*yyprod(:).*gg)'*da_stacked*dz;
% d = max(ddprod.*(kk-aa),0);
% Debt = (d(:).*gg)'*da_stacked*dz;
% DY=Debt/GDP;

% CREATE PLOTS COMPARING RESULTS 
load tax_iteration_results_lf.mat;
load tax_iteration_results_opt.mat;

% Figures with optimal credit subsidy as deviation from no credit subsidy
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

% Figures with both optimal credit subsidy and no credit subsidy
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,taulvec_opt,'-b');
xlabel('Time');
title(sprintf('Optimal credit subsidy schedule, half life = %5.3f',halflife_opt));
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_taxschedule%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,wage_vec_opt,'-b',timevec,wage_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Wage');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_wage%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,profits_total_vec_opt,'-b',timevec,profits_total_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Total profits');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_profitstot%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,profit_share_vec_opt,'-b',timevec,profit_share_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Profits/(profits + labor income)');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_profitshare%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_entrepreneurs_opt,'-b',timevec,frac_entrepreneurs_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Fraction entrepreneurs');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracentrepreneurs%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_constrained_vec_opt,'-b',timevec,frac_constrained_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Fraction entrepreneurs constrained');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracconstrained%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_output_constrained_vec_opt,'-b',timevec,frac_output_constrained_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Fraction output from constrained entrepreneurs');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracoutputconstrained%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,labor_vec_opt,'-b',timevec,labor_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Aggregate labor supply');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_labor%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,capital_vec_opt,'-b',timevec,capital_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Aggregate capital use');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_capital%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,output_vec_opt,'-b',timevec,output_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('GDP');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_gdp%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,tfp_vec_opt,'-b',timevec,tfp_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('TFP');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_tfp%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,frac_output_of_unconstrained_vec_opt,'-b',timevec,frac_output_of_unconstrained_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('GDP/Potential GDP');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_fracgdpofpotential%d',iteration)]);
delete(H);
 
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,wealth_tot_vec_opt,'-b',timevec,wealth_tot_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Aggregate wealth');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_aggwealth%d',iteration)]);
delete(H);

H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,utility_vec_opt,'-b',timevec,utility_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Period utility');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_periodutility%d',iteration)]);
delete(H);
   
H=figure;
set(gcf,'Visible','off');
xlim([0 T]);
plot(timevec,wage_vec_opt.*labor_vec_opt,'-b',timevec,wage_vec_lf.*labor_vec_lf,'-r');
legend('Optimal credit subsidy','No credit subsidy');
xlabel('Time');
title('Worker period consumption');
set(gca,'FontSize',14);
print(print_graph_type,[path_graphs,'/',sprintf('fig_periodconsumption%d',iteration)]);
delete(H);

end

save plotting_vectors.mat timevec labor_vec_opt labor_vec_lf wealth_tot_vec_opt wealth_tot_vec_lf output_vec_opt output_vec_lf wage_vec_opt wage_vec_lf ...
tfp_vec_opt tfp_vec_lf utility_vec_opt utility_vec_lf taulvec_opt profits_total_vec_opt profits_total_vec_lf profit_share_vec_opt profit_share_vec_lf ...
frac_entrepreneurs_opt frac_entrepreneurs_lf frac_constrained_vec_opt frac_constrained_vec_lf frac_output_constrained_vec_opt frac_output_constrained_vec_lf ...
capital_vec_opt capital_vec_lf frac_output_of_unconstrained_vec_opt frac_output_of_unconstrained_vec_lf DY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE LATEX FILE        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVE PLOTS TO LATEX FILE
if print_optpath==1
    
load tax_iteration_results_lf.mat;
load tax_iteration_results_opt.mat;
load tax_iteration_results_optss.mat;

copyfile('output.tex','texfile.tex');
fid=fopen('texfile.tex','a');

% Table comparing welfare of various credit subsidy schedules

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
fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | c | c | c } \n \\hline \n');
fprintf(fid,' & No cred sub & Optimal policy & Optimal SS cred sub & Optimal Flat cred sub \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'$\\varsigma_0$ & %7.5f & %7.5f & %7.5f & %7.5f \\\\ \n',taulvec_lf(1),taulvec_opt(1),taulvec_optss(1),opttaul_flat);
fprintf(fid,'$\\bar{\\varsigma}$ & %7.5f & %7.5f & %7.5f & %7.5f \\\\ \n',taulvec_lf(end),taulvec_opt(end),taulvec_optss(end),opttaul_flat);
fprintf(fid,'Half life & - & %7.5f & - & - \\\\ \n',halflife_opt);
fprintf(fid,'Welfare (weighted) & %7.5f & %7.5f & %7.5f & %7.5f \\\\ \n',welfare_lf,welfare_opt,welfare_optss,welfare_flat);
fprintf(fid,'Welfare workers & %7.5f & %7.5f & %7.5f & %7.5f \\\\ \n',welfare_workers_lf,welfare_workers_opt,welfare_workers_optss,welfare_workers_flat);
fprintf(fid,'Welfare entrepreneurs & %7.5f & %7.5f & %7.5f & %7.5f \\\\ \n',welfare_entrepreneurs_lf,welfare_entrepreneurs_opt,welfare_entrepreneurs_optss,welfare_entrepreneurs_flat);
fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n');
fprintf(fid,'\\end{table} \n \n');

fprintf(fid,'\\begin{table}[h!]\n');
fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | c } \n \\hline \n');
fprintf(fid,' & Constant $\\varsigma_0$ & Constant $\\bar{\\varsigma}$ \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'$\\varsigma_0$ & %7.5f & %7.5f \\\\ \n',taulvec_t0(1),taulvec_opt(end));
fprintf(fid,'$\\bar{\\varsigma}$ & %7.5f & %7.5f \\\\ \n',taulvec_t0(end),taulvec_opt(end));
fprintf(fid,'Half life & - & - \\\\ \n');
fprintf(fid,'Welfare (weighted) & %7.5f & %7.5f \\\\ \n',welfare_t0,welfare_bar);
fprintf(fid,'Welfare workers & %7.5f & %7.5f \\\\ \n',welfare_workers_t0,welfare_workers_bar);
fprintf(fid,'Welfare entrepreneurs & %7.5f & %7.5f \\\\ \n',welfare_entrepreneurs_t0,welfare_entrepreneurs_bar);
fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n');
fprintf(fid,'\\end{table} \n \n');

%fprintf(fid,'Optimal flat credit subsidy rate = %5.3f \n',opttaul_flat);

fprintf(fid,'\\clearpage \n \n');

rho_ave = (pop_share/discrate_workers+(1-pop_share)/rho)^(-1);

fprintf(fid,'\\begin{table}[h!]\n');
fprintf(fid,'\\begin{center} \n \\begin{tabular}{l | c | c | c } \n \\hline \n');
fprintf(fid,' Experiment & Total welfare & Worker welfare & Entrepreneur welfare \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'Optimal policy & %7.5f & %7.5f & %7.5f \\\\ \n',exp(rho_ave*(welfare_opt-welfare_lf))-1,exp(discrate_workers*(welfare_workers_opt-welfare_workers_lf))-1,exp(rho*(welfare_entrepreneurs_opt-welfare_entrepreneurs_lf))-1);
fprintf(fid,'Optimal flat cred sub & %7.5f & %7.5f & %7.5f \\\\ \n',exp(rho_ave*(welfare_flat-welfare_lf))-1,exp(discrate_workers*(welfare_workers_flat-welfare_workers_lf))-1,exp(rho*(welfare_entrepreneurs_flat-welfare_entrepreneurs_lf))-1);
fprintf(fid,'Constant $\\varsigma_0$ & %7.5f & %7.5f & %7.5f \\\\ \n',exp(rho_ave*(welfare_t0-welfare_lf))-1,exp(discrate_workers*(welfare_workers_t0-welfare_workers_lf))-1,exp(rho*(welfare_entrepreneurs_t0-welfare_entrepreneurs_lf))-1);
fprintf(fid,'Constant $\\bar{\\varsigma}$ & %7.5f & %7.5f & %7.5f \\\\ \n',exp(rho_ave*(welfare_bar-welfare_lf))-1,exp(discrate_workers*(welfare_workers_bar-welfare_workers_lf))-1,exp(rho*(welfare_entrepreneurs_bar-welfare_entrepreneurs_lf))-1);
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
fprintf(fid,'\\item Credit subsidy schedule: $ \\varsigma_k(t) = \\bar{\\varsigma}_k + e^{-\\gamma t} (\\varsigma_{k,0}-\\bar{\\varsigma}_k) $ \n');
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
fprintf(fid,'Range of initial credit subsidy rate tested & $\\varsigma_0$ & [%5.3f,%5.3f] \\\\ \n',tau0_vec(1),tau0_vec(end));
fprintf(fid,'Range of final credit subsidy rate tested & $\\bar{\\varsigma}$ & [%5.3f,%5.3f] \\\\ \n',taubar_vec(1),taubar_vec(end));
fprintf(fid,'Contraction of initial distribution & $factor$ & %5.3f \\\\ \n',factor);
fprintf(fid,'\\hline \n \\end{tabular} \n \\end{center} \n \n');

fprintf(fid,'\\clearpage \n \n');

% Section: figures
fprintf(fid,'\\section{Figures} \n');

end

% Print Optimal Credit Subsidy Rate Figures
if print_optss==1

fprintf(fid,'Optimal steady state credit subsidy rate = %5.3f \n',opttaul_ss);
% Welfare as function of credit subsidy
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_welfaress%d',iteration)]);
fprintf(fid,'\\caption{}\n');
fprintf(fid,'\\end{figure}\n\n'); 
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

% PRINT Welfare by credit subsidy schedule
if print_findpath==1
    % Print 3D plots of welfare by credit subsidy schedule parameters
    for ii=1:length(halflives_vec)
        % Relative welfare: workers
        fprintf(fid,'\\begin{figure}[h!]\n');
        fprintf(fid,'\\centering\n');
        fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_welfarepathtax_workers_hl%d_%d',halflives_vec(ii),iteration)]);
        fprintf(fid,'\\caption{}\n');
        fprintf(fid,'\\end{figure}\n\n'); 
        % Relative welfare: entrepreneurs
        fprintf(fid,'\\begin{figure}[h!]\n');
        fprintf(fid,'\\centering\n');
        fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_welfarepathtax_entrepreneurs_hl%d_%d',halflives_vec(ii),iteration)]);
        fprintf(fid,'\\caption{}\n');
        fprintf(fid,'\\end{figure}\n\n'); 
        % Total welfare: total
        fprintf(fid,'\\begin{figure}[h!]\n');
        fprintf(fid,'\\centering\n');
        fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_welfarepathtax_hl_total%d_%d',halflives_vec(ii),iteration)]);
        fprintf(fid,'\\caption{}\n');
        fprintf(fid,'\\end{figure}\n\n'); 
    end
end

fprintf(fid,'\\clearpage \n \n');

% Figures with optimal credit subsidy as deviation from no credit subsidy
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_deviationoptlf%d',iteration)]);
fprintf(fid,'\\caption{Proportional deviations of optimal credit subsidy equilibrium from the laissez-faire equilibrium}\n');
fprintf(fid,'\\end{figure}\n\n'); 

% Figures with both no credit subsidy and optimal credit subsidy
fprintf(fid,'\\begin{figure}[h!]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',[path_graphs_rel,'/',sprintf('fig_taxschedule%d',iteration)]);
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
copyfile('main_Credit_Supply_0.m',path_docs);

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
