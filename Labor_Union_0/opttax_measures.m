% PRINTS VARIOUS MEASURES FROM TRANSITION EQUILIBRIUM
% Inputs: w_t, taulvec
% opttax_noc_parameters.m must be run first

if preftype==1
    transfervec = taulvec.*(1-taulvec).^chi0.*w_t.^(1+chi0);
    lsvec = ((1-taulvec).*w_t).^chi0;
    consumption_workers = w_t.*lsvec;
    if gam==1
        utility_vec = log((consumption_workers)-(lsvec).^(1+1/chi0)/(1+1/chi0));
        %utility_vec2 = log((1+chi0)^(-1)*((1-taulvec).*w_t).^(1+chi0)+transfervec);
    else
        utility_vec = (1-gam)^(-1)*((consumption_workers)-(lsvec).^(1+1/chi0)/(1+1/chi0)).^(1-gam); 
        %utility_vec2 = (1-gam)^(-1)*((1+chi0)^(-1)*((1-taulvec).*w_t).^(1+chi0)+transfervec).^(1-gam);
    end
elseif preftype==2
    lsvec=(w_t.^(1-gam).*(1-taulvec)).^(1/(gam+1/chi0));
    consumption_workers = w_t.*lsvec;
    if gam==1
        utility_vec = (log(consumption_workers)-(lsvec).^(1+1/chi0)./(1+1/chi0));
    else
        utility_vec = ((consumption_workers)^(1-gam)/(1-gam)-(lsvec).^(1+1/chi0)./(1+1/chi0));
    end
end

consumption_workers_vec = consumption_workers;
profits_total_vec = zeros(N,1);
profit_share_vec = zeros(N,1);
frac_active_vec= zeros(N,1);
frac_constrained_vec = zeros(N,1);
frac_output_constrained_vec = zeros(N,1);
frac_output_of_unconstrained_vec = zeros(N,1);
wealth_tot_vec = zeros(N,1);
capital_vec = zeros(N,1);
labor_vec = zeros(N,1);
output_vec = zeros(N,1);
tfp_vec = zeros(N,1);
frac_entrepreneurs = (1-pop_share)*ones(N,1);
debt_vec = zeros(N,1);
dy = zeros(N,1);

for ii=1:N
    w=w_t(ii);
    taul=taulvec(ii);

    kk=min(la*aa,(zz*Aprod).^(1/(1-al-be))*(al/(rstar+d))^((1-be)/(1-al-be))*...
        (be/w)^(be/(1-al-be))+fk); %Capital demand
    nn=(be*zz*Aprod/w).^(1/(1-be)).*max(kk-fk,0).^(al/(1-be))+fn; %Labor demand
    yyprod=zz.*Aprod.*max(kk-fk,0).^al.*max(nn-fn,0).^be; %Output
    Pi=yyprod-(rstar+d).*kk-w*nn; %Profits
    ddprod = Pi>0; %Indicator for whether entrepreneur is active
    M=ddprod.*Pi; %Income of entrepreneurs
    
    kk_unconstrained = (zz*Aprod).^(1/(1-al-be))*(al/(rstar+d))^((1-be)/(1-al-be))*...
    (be/w)^(be/(1-al-be))+fk;
    nn_unconstrained=(be*zz*Aprod/w).^(1/(1-be)).*max(kk_unconstrained-fk,0).^(al/(1-be))+fn;
    yyprod_unconstrained=zz.*Aprod.*max(kk_unconstrained-fk,0).^al.*max(nn_unconstrained-fn,0).^be;
    
    % Profits total
    profits_total_vec(ii) = (1-pop_share)*((M(:).*gg{ii})'*da_stacked*dz);
    % Profit share of income
    profit_share_vec(ii) = ((M(:).*gg{ii})'*da_stacked*dz)/((M(:).*gg{ii})'*da_stacked*dz+w*ls_it(ii,it));
    % Fraction of entrepreneurs active
    frac_active_vec(ii) = (ddprod(:).*gg{ii})'*da_stacked*dz;
    % Fraction of entrepreneurs constrained
    frac_constrained_vec(ii) = (ddconstrained(:).*gg{ii})'*da_stacked*dz;
    % Fraction of output produced by constrained entrepreneurs
    frac_output_constrained_vec(ii) = ((yyprod(:).*ddprod(:).*ddconstrained(:).*gg{ii})'*da_stacked*dz)/((yyprod(:).*ddprod(:).*gg{ii})'*da_stacked*dz);
    % Actual output as a fraction of total output if no-one is constrained
    frac_output_of_unconstrained_vec(ii) = ((ddprod(:).*yyprod(:).*gg{ii})'*da_stacked*dz)/((yyprod_unconstrained(:).*gg{ii})'*da_stacked*dz);
    % Total wealth
    wealth_tot_vec(ii) = (1-pop_share)*(as.*gg{ii})'*da_stacked*dz;
    
    % TFP
    capital_vec(ii) = (1-pop_share)*(ddprod(:).*kk(:).*gg{ii})'*da_stacked*dz;
    labor_vec(ii) = (1-pop_share)*(ddprod(:).*nn(:).*gg{ii})'*da_stacked*dz;
    output_vec(ii) = (1-pop_share)*(ddprod(:).*yyprod(:).*gg{ii})'*da_stacked*dz;
    tfp_vec(ii) = output_vec(ii)/(capital_vec(ii)^al*labor_vec(ii)^be);
    debt_vec(ii) = (1-pop_share)*(ddprod(:).*max(kk(:)-aa(:),0).*gg{ii})'*da_stacked*dz;
    dy(ii) = debt_vec(ii)/output_vec(ii);
end

% Frisch elasticity
if preftype==1
    uarg = ((consumption_workers)-(lsvec).^(1+1/chi0)/(1+1/chi0));
    uc = uarg.^(-gam);
    ul = -uarg.^(-gam).*(lsvec).^(1/chi0);
    ucc = -gam*uarg.^(-gam-1);
    ull = uarg.^(-gam-1).*(lsvec).^(2/chi0)-(1/chi0)*uarg.^(-gam).*(lsvec).^(1/chi0-1);
    ucl = gam*uarg.^(-gam-1).*(lsvec).^(1/chi0);
    frisch_el = -((-ul./(lsvec)).*ucc)./(ucl.^2-ull.*ucc);
    frisch_el = -((-ul./(lsvec))./(ucl.^2./ucc-ull)); % without tax = consant = chi0/(2*(1+chi0)*chi0-1)
end
