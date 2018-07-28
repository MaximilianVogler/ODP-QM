% disp('Iteration    Wage    VF steps    DIST steps    Excess supply');
flag=0;

%INITIALIZE PRICES
wmin=wmin0; wmax = wmax0; w=w0;
transfer=0;

%START ITERATION OVER WAGE HERE
for it=1:maxit_price
% fprintf('    %d       ',it);
% fprintf('%7.5f      ',w);

%SOLVE STATIC PROBLEM OF CONSUMPTION AND LABOR
kk=min(la*aa,(zz*Aprod).^(1/(1-al-be))*(al/(rstar+d))^((1-be)/(1-al-be))*...
    (be/w)^(be/(1-al-be))+fk); %Capital demand
ddconstrained = (kk==la*aa);
nn=(be*zz*Aprod/w).^(1/(1-be)).*max(kk-fk,0).^(al/(1-be))+fn; %Labor demand
yyprod=zz.*Aprod.*max(kk-fk,0).^al.*max(nn-fn,0).^be; %Output
Pi=yyprod-(rstar+d).*kk-w*nn; %Profits
ddprod = Pi>0;
M=ddprod.*Pi; %Income of entrepreneurs

kk_unconstrained = (zz*Aprod).^(1/(1-al-be))*(al/(rstar+d))^((1-be)/(1-al-be))*...
    (be/w)^(be/(1-al-be))+fk;
nn_unconstrained=(be*zz*Aprod/w).^(1/(1-be)).*max(kk_unconstrained-fk,0).^(al/(1-be))+fn;
yyprod_unconstrained=zz.*Aprod.*max(kk_unconstrained-fk,0).^al.*max(nn_unconstrained-fn,0).^be;

%Guess value function
if it==1
    if gam==1
        v0=log(M+rstar*aa)/(rho+the);
    else
        v0 = (M+rstar*aa).^(1-gam)/(1-gam)/(rho+the);
    end
    v=v0;
end

% ITERATION OVER VALUE FUNCTION
for n=maxit_vf:-1:1
    %v is the old value function (n+1) from previous step
    V = v; 
    % forward difference
    Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./daaf(1:I-1,:);
    Vaf(I,:) = 0; %will never be used
    % backward difference
    Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))./daab(2:I,:);
    Vab(1,:) = (M(1,:) + rstar.*amin).^(-gam); %state constraint boundary condition
    
    I_concave = Vab > Vaf; %indicator whether value function is concave (problems arise if this is not the case)
    
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
    AA = -A + (1/Delta + rho + the)*speye(I*J);
    
   
    %SOLVE HAMILTON-JACOBI-BELMAN SYSTEM
    u_stacked = reshape(u,I*J,1); %column vector, first z1 for all a, then z2...
    V_stacked = reshape(V,I*J,1);
    
    b = u_stacked + V_stacked/Delta;

    V_stacked = AA\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = reshape(V_stacked,I,J);
    
    Vchange = V - v;
    v = V;

    dist = max(max(abs(Vchange)));
    if dist<crit_vf
        break;
    elseif n==1
        disp('No convergence of value function!');
        flag=1;
        break;
    end
end

%fprintf('%d          ',maxit_vf-n+1);


%KOLMOGOROV FORWARD EQUATION
%ADJUST TO ENSURE MASS-PRESERVATION WITH NON-UNIFORM GRID
grid_diag = spdiags(da_stacked,0,I*J,I*J);
Anew = grid_diag*A*grid_diag^(-1);
AT = Anew';

% clear gg;
% gg{1}=gg_sparse;
% %gg{1}=gg0;
% for m=1:maxit_dist
%     amean=(as.*gg{m})'*da_stacked*dz;
%     zmean=(zs.*gg{m})'*da_stacked*dz;
%     Psi=zeros(I,J);
%     [aux,zmeanind]=min(abs(z-zmean));
%     ind1=find(a<amean,1,'last');
%     weight1=(a(ind1+1)-amean)/da_tilde(ind1+1);
%     Psi(ind1,zmeanind)=weight1*(dz*da_tilde(ind1+1))^(-1);
%     Psi(ind1+1,zmeanind)=(1-weight1)*(dz*da_tilde(ind1+1))^(-1);   
%     Psi=sparse(reshape(Psi,I*J,1));
%     Psi=Psi/(Psi(:)'*da_stacked*dz);
%     %Implicit updating of gg
%     gg{m+1}= (speye(J*I) - AT*Delta+the*Delta*speye(J*I))\(gg{m}+the*Delta*Psi);
%     gdiff=max(abs(gg{m+1}-gg{m}));
%     %Check for convergence
%     if gdiff<crit_dist
%         break;
%     elseif m==maxit_dist
%         disp('No convergence of distribution!');
%         flag=1;
%         break;
%     end
% end
% gg_sparse=gg{m};
% gg=full(gg{m});

i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
AT(i_fix,:) = row;
gg = sparse(AT)\sparse(b);
gg = full(gg);


%NORMALIZE AND PREPARE MARGINALS
g_sum = gg'*da_stacked*dz; %Normalize density
gg = gg./g_sum;
g = reshape(gg,I,J);
g_z = da_tilde'*g;
g_a = sum(g.*dz,2);


%UPDATE WAGE AND CHECK CONVERGENCE
if preftype==1
    ls = pop_share*((1-taul)*w)^chi0; %Supply of labor (total)
elseif preftype==2
    ls = pop_share*(w^(1-gam)*(1-taul))^(1/(gam+1/chi0));
end
ld = (1-pop_share)*(ddprod(:).*nn(:).*gg)'*da_stacked*dz; %Demand of labor
ES=ls-ld; %Excess supply

% Update wage
if ES>crit_price
    wmax=w;
    w = 0.5*(wmax+wmin);
elseif ES<-crit_price
    wmin = w;
    w = 0.5*(wmax+wmin);
elseif abs(ES)<crit_price
    break;
end

%fprintf('%d    %7.5f    %7.5f    %7.5f \n',it,w,ES,transferdiff);

%Break if either VF or distribution did not converge
if flag==1, break; end;
%Break if exceed max iteration on wage
if it==maxit_price
    disp('No wage convergence!');
end


end

%LUMP-SUM TRANSFER AND WELFARE
%transfer = taul*((1-taul)*w)^chi0;
ls=ls/pop_share; % Per capita labor
welfare_entrepreneurs = (v(:).*gg)'*da_stacked*dz;  % Per capita entrepreneur welfare
consumption_workers = w*ls; %Per capita worker consumption level (=income)
%transfer = taul*w*ls;
if preftype==1
    %labordisutility=((1-taul)*w)^(1+1/chi0)/(1+1/chi0);
    %transfer = taul*(1-taul)^chi0*w^(1+chi0);
    if gam==1
        welfare_workers = (discrate_workers)^(-1)*log((consumption_workers)-(ls).^(1+1/chi0)/(1+1/chi0));  % Per capita welfare of worker
        %welfare_workers2 = L*(discrate_workers)^(-1)*log((1+chi0)^(-1)*((1-taul)*w)^(1+chi0)+transfer);
    else
        welfare_workers = (discrate_workers)^(-1)*(1-gam)^(-1)*((consumption_workers)-(ls).^(1+1/chi0)/(1+1/chi0)).^(1-gam);
        %welfare_workers2 = L*((1-gam)*(discrate_workers))^(-1)*((1+chi0)^(-1)*((1-taul)*w)^(1+chi0)+transfer)^(1-gam);
    end
elseif preftype==2
    if gam==1
        welfare_workers = (discrate_workers)^(-1)*(log(consumption_workers)-(ls)^(1+1/chi0)/(1+1/chi0));
    else
        welfare_workers = (discrate_workers)^(-1)*((consumption_workers)^(1-gam)/(1-gam)-(ls)^(1+1/chi0)/(1+1/chi0));
    end
end

%gg=full(gg_sparse);

% VARIOUS MEASURES
% Profit share of income
profits_tot = (1-pop_share)*(M(:).*gg)'*da_stacked*dz;
wages_tot = pop_share*w*ls;
profit_share = profits_tot/(profits_tot+wages_tot);
% Fraction of entrepreneurs active
frac_active = (ddprod(:).*gg)'*da_stacked*dz;
% Fraction of entrepreneurs constrained
frac_constrained = (ddconstrained(:).*gg)'*da_stacked*dz;
% Fraction of output produced by constrained entrepreneurs
frac_output_constrained = ((yyprod(:).*ddprod(:).*ddconstrained(:).*gg)'*da_stacked*dz)/((yyprod(:).*ddprod(:).*gg)'*da_stacked*dz);
% Actual output as a fraction of total output if none is constrained
frac_output_of_unconstrained = ((ddprod(:).*yyprod(:).*gg)'*da_stacked*dz)/((yyprod_unconstrained(:).*gg)'*da_stacked*dz);
% Total wealth
wealth_tot = (1-pop_share)*(as.*gg)'*da_stacked*dz;

% Frisch elasticity
if preftype==1
    uarg = consumption_workers-(ls).^(1+1/chi0)/(1+1/chi0);
    uc = uarg.^(-gam);
    ul = -uarg.^(-gam).*(ls).^(1/chi0);
    ucc = -gam*uarg.^(-gam-1);
    ull = uarg.^(-gam-1).*(ls).^(2/chi0)-(1/chi0)*uarg.^(-gam).*(ls).^(1/chi0-1);
    ucl = gam*uarg.^(-gam-1).*(ls).^(1/chi0);
    frisch_el_ss = -((-ul./(ls)).*ucc)./(ucl.^2-ull.*ucc);
end
