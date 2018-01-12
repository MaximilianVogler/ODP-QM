% Same as v4, but prepared for parfor loop over taxes
% Put inside a function wrapping because this is needed for parfor

function [welfare_mat,gg,w_t,pi_t,flag,various] = opttax_noc_transition(gg_guess,gg_initial,taulvec)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL VALUE FUNCTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run opttax_noc_parameters;
%load test.mat;
gg_sparse = sparse(gg_guess); %Guess distribution from steady state without tax
taul=taulvec(N);
run opttax_noc_steadystate_no_LU.m;
v_st = v; %Used as final value function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSITION PATH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run opttax_noc_parameters;
wnew = w*ones(N,1); %Guess wage sequence
transfervec=zeros(N,1); %Guess lump-sum transfers
    
%Initialize matrices for storing transition path
v = zeros(I,J,N);
gg = cell(N,1);
g_it = zeros(I,J,N);
g_a_it = zeros(I,N);
g_z_it = zeros(J,N);
lab_constr = zeros(N,1);
v(:,:,N)= v_st;
A_t=cell(N+1,1);
nnn = zeros(I,J,N);
dddprod = zeros(I,J,N);
dddconstrained = zeros(I,J,N);
w_it=zeros(N,maxit_price_seq);
pi_t=zeros(N,1);
ES_it=zeros(N,maxit_price_seq); %Saves iterations on the path of excess supply as columns
ls_it=zeros(N,maxit_price_seq);
ld_it=zeros(N,maxit_price_seq);
%s_it=zeros(I,J,N);

%ITERATION OVER WAGE SEQUENCE STARTS HERE
for it=1:maxit_price_seq

w_t = wnew; %Update wage sequence
w_it(:,it)=w_t;

%SOLVE FOR PATH OF VALUE FUNCTIONS
V = v_st;
for n=N:-1:1
    v(:,:,n)=V;
    w=w_t(n);
    %taul=taulvec(n);
    
    %SOLVE STATIC PROBLEM OF CONSUMPTION AND LABOR
    kk=min(la*aa,(zz*Aprod).^(1/(1-al-be))*(al/(rstar+d))^((1-be)/(1-al-be))*...
        (be/w)^(be/(1-al-be))+fk); %Capital demand
    ddconstrained = (kk==la*aa);
    nn=(be*zz*Aprod/w).^(1/(1-be)).*max(kk-fk,0).^(al/(1-be))+fn; %Labor demand
    yyprod=zz.*Aprod.*max(kk-fk,0).^al.*max(nn-fn,0).^be; %Output
    Pi=yyprod-(rstar+d).*kk-w*nn; %Profits
    ddprod = Pi>0; %Indicator for whether entrepreneur is active
    M=ddprod.*Pi; %Income of entrepreneurs
    
    dddconstrained(:,:,n) = ddconstrained;
    nnn(:,:,n) = nn; %Save labor choice for later use
    dddprod(:,:,n) = ddprod; 

    % forward difference
    Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./daaf(1:I-1,:);
    Vaf(I,:) = 0; %will never be used
    % backward difference
    Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))./daab(2:I,:);
    %Vab(1,:) = (M(1,:) + rstar.*amin).^(-gam); %state constraint boundary condition
    Vab(1,:) = (M(1,:) + rstar.*amin).^(-gam); %state constraint boundary condition
    
    %I_concave = Vab > Vaf; %indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf = Vaf.^(-1/gam);
    sf = M + rstar.*aa - cf;
    %consumption and savings with backward difference
    cb = Vab.^(-1/gam);
    sb = M + rstar.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = M + rstar.*aa;
    Va0 = (c0).^(-gam);
    
    %upwind scheme   
    If = sf > 0; %positive drift --> forward difference
    Ib = sb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    
    Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term
    
    c = Va_Upwind.^(-1/gam);
    if gam==1
        u = log(c);
    else
        u = (c).^(1-gam)/(1-gam);
    end
    s = M + rstar.*aa - c;
    
    %CONSTRUCT MATRIX A
    X = -min(sb,0)./daab;
    Y = -max(sf,0)./daaf + min(sb,0)./daab;
    Z = max(sf,0)./daaf;
    
    updiag=0; %This is needed because of the peculiarity of spdiags.
    for j=1:J
        updiag=[updiag;Z(1:I-1,j);0];
    end
    
    centdiag=reshape(Y,I*J,1);
    
    lowdiag=X(2:I,1);
    for j=2:J
        lowdiag=[lowdiag;0;X(2:I,j)];
    end
    
    B=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
    
    A = B + Bswitch;
    A_t{n} = A;
    AA = -A + (1/Deltavec(n) + rho)*speye(I*J);
    
    %SOLVE HAMILTON-JACOBI-BELMAN SYSTEM
    u_stacked = reshape(u,I*J,1); %column vector, first z1 for all a, then z2...
    V_stacked = reshape(V,I*J,1);
    
    b = u_stacked + V_stacked/Deltavec(n);

    V_stacked = AA\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = reshape(V_stacked,I,J);
end

gg{1}=sparse(gg_initial);
lab_constr = zeros(N,1);
for m=1:N
    Anew = grid_diag*A_t{m}*grid_diag^(-1);
    AT = Anew';
    amean=(as.*gg{m})'*da_stacked*dz;
    zmean=(zs.*gg{m})'*da_stacked*dz;
    Psi=zeros(I,J);
    [aux,zmeanind]=min(abs(z-zmean));
    ind1=find(a<amean,1,'last');
    weight1=(a(ind1+1)-amean)/da_tilde(ind1+1);
    Psi(ind1,zmeanind)=weight1*(dz*da_tilde(ind1+1))^(-1);
    Psi(ind1+1,zmeanind)=(1-weight1)*(dz*da_tilde(ind1+1))^(-1);   
    Psi=sparse(reshape(Psi,I*J,1));
    Psi=Psi/(Psi(:)'*da_stacked*dz);
    %Implicit updating of gg
    gg{m+1} = (speye(J*I) - AT*Deltavec(m+1)+the*Deltavec(m+1)*speye(J*I))\(gg{m}+the*Deltavec(m+1)*Psi);
    gg{m+1} = gg{m+1}/(gg{m+1}'*da_stacked*dz); %Not really needed if gg{1} has total mass 1, since iteration is mass-preserving
    %Excess labor supply
    ddprod = reshape(dddprod(:,:,m),I*J,1);
    nns = reshape(nnn(:,:,m),I*J,1);
    if preftype==1
        ls = pop_share*((1-taulvec(m))*w_t(m)).^chi0; %Supply of labor (total)
    elseif preftype==2
        ls = pop_share*(w_t(m)^(1-gam)*(1-taulvec(m)))^(1/(gam+1/chi0));
    end
    ld = (1-pop_share)*(ddprod(:).*nns(:).*gg{m})'*da_stacked*dz; %Demand of labor
    ES_it(m,it)=ls-ld; %Excess supply
    ls_it(m,it)=ls;
    ld_it(m,it)=ld;
end

% Check convergence of pi
for n=1:N
    ddddconstrained=reshape(dddconstrained(:,:,n),I*J,1);
    nnnn=reshape(nnn(:,:,n),I*J,1);
    ddprod = reshape(dddprod(:,:,n),I*J,1);
    lab_constr(n) = (1-pop_share)*(ddprod(:).*(1-ddddconstrained(:)).*nnnn(:).*gg{n})'*da_stacked*dz;
    pi_t(n) = lab_constr(n)/ld_it(n,it);
end

% UPDATE TRANSFER
% if preftype==2
%     transfervec_new = taulvec.*w_t.*ls_it(:,it);
%     transferdiff = transfervec_new-transfervec;
%     [~,indaux] = max(abs(transferdiff));
%     maxtransferdiff = transferdiff(indaux);
%     transfervec=transfervec_new;
% end

%UPDATE WAGE
if it<10
    adj = 0.2; %0.2;
elseif it<20
    adj = 0.1;
elseif it<30
    adj = 0.1;
elseif it<40
    adj = 0.05;
else
    adj = 0.01;
end
wnew = w_t.*(1-min(adj*ES_it(:,it),0.9)); % Adjust based on ES (ES>0, then w should decrease)

%CHECK CONVERGENCE
maxwagediff = max(abs(wnew-w_t));
[~,indaux] = max(abs(ES_it(:,it)));
maxESdiff=ES_it(indaux,it);
% Break if maximum ES low enough
if abs(maxESdiff)<crit_price_seq
    break;
end
% Break if maximum number of iterations reached
if it==maxit_price_seq
    disp('No convergence of excess labor supply'); %Print error if no convergence
    flag=2;
end


fprintf('    %d           %8.6f       %8.6f \n',it,maxwagediff,maxESdiff);


end


% TRANSITION STATISTICS
%Distributions
for n=1:N
    g = reshape(gg{n},I,J);
    g_it(:,:,n) = g;
    g_z_it(:,n) = da_tilde'*g;
    g_a_it(:,n) = sum(g.*dz,2);
end


% WELFARE
%Lump-sum transfers and period utility
%transfervec = L*taulvec.*(1-taulvec).^chi0.*w_t.^(1+chi0);
lsvec=ls_it(:,it)/pop_share; %Per capita labor supply
consumption_workers = w_t.*lsvec; %Total consumption of workers
%Find vector of per-worker utility levels by period
if preftype==1
    if gam==1
        utility_vec = log(consumption_workers-(lsvec).^(1+1/chi0)/(1+1/chi0)); %log((1+chi0)^(-1)*((1-taulvec).*w_t).^(1+chi0)+transfervec);    % Welfare per capita
    else
        utility_vec = (1-gam)^(-1)*(consumption_workers-(lsvec).^(1+1/chi0)/(1+1/chi0)).^(1-gam); %(1-gam)^(-1)*((1+chi0)^(-1)*((1-taulvec).*w_t).^(1+chi0)+transfervec).^(1-gam);
    end
elseif preftype==2
    if gam==1
        utility_vec = (log(consumption_workers)-(lsvec).^(1+1/chi0)./(1+1/chi0));
    else
        utility_vec = ((consumption_workers)^(1-gam)/(1-gam)-(lsvec).^(1+1/chi0)./(1+1/chi0));
    end
end

% Two trapezoidal approximations of welfare integral over time (the second
% uses an exact integration of the exponential function within periods
disc_vec = exp(-discrate_workers*timevec');
Deltavec2 = Deltavec(2:end-1);
dt = zeros(N,1); 
dtf = zeros(N,1); 
dtb = zeros(N,1);
dtf(1:N-1) = timevec(2:N) - timevec(1:N-1);
dtb(2:N) = timevec(2:N) - timevec(1:N-1);
dt = 0.5*(dtb + dtf);
dt(1) = 0.5*dtf(1); dt(N) = 0.5*dtb(N);

%Total worker welfare
welfare_workers_trapezoidal1 = ( Deltavec(1)*utility_vec(1)/2 ...
    +(Deltavec2(2:end)'+Deltavec2(1:end-1)')/2.*exp(-discrate_workers*timevec(2:end-1)')*utility_vec(2:end-1) ...
    +Deltavec2(end)*exp(-discrate_workers*T)*utility_vec(end)/2 );
welfare_workers_trapezoidal2 = discrate_workers^(-1)*( utility_vec(1)/2*(disc_vec(1)-disc_vec(2)) +...
    (disc_vec(1:end-2)-disc_vec(3:end))*utility_vec(2:end-1)/2 + utility_vec(end)/2*(disc_vec(end-1)-disc_vec(end)) );
welfare_workers_sum = sum(dt.*exp(-discrate_workers*timevec).*utility_vec) + exp(-discrate_workers*timevec(end))*utility_vec(end)/discrate_workers; 

% Entrepreneurs' welfare
welfare_entrepreneurs = (reshape(v(:,:,1),I*J,1).*gg{1})'*da_stacked*dz; %Period 0 welfare with zero weight on future generations

% Matrix of welfare
welfare_mat = [welfare_workers_trapezoidal1 welfare_workers_sum welfare_entrepreneurs];

% Output
various = ES_it; %ls_it;

end