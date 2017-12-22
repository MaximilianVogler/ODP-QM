%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMON PARAMETERS AND INITIALIZATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREFERENCES AND PRODUCTION
paretoweight_workers = 4/5; %Pareto weight on workers
%paretoweight_ent =1/3;
pop_share = 2/3; %Population share of workers
pop_mass = 1; %Mass of entire population, i.e. sum of entrepeneurs and workers
L = pop_share*pop_mass; %Mass of workers (relative to entrepreneurs)
preftype = 2; % Type of utility function for workers; 1=GHH, 2=balanced growth
poissonshocks = 1; %Indicator for whether there are Poisson shocks to productivity
poissonarrivalrate = 0.1; %Arrival rate of Poisson shocks
gam = 1; %CRRA utility with parameter gamma
chi0 = 1; %Labor disutility parameter (~ Frisch elasticity)
d = 0; %depreciation rate
rho = 0.05; %0.1; %0.05; %discount rate for entrepreneurs
the = 0; %death rate
fk = 0; %2; %fixed capital cost
fn = 0; %2; %fixed labor cost
la = 2; %financial constrain parameter
eta = 0; %Effect of productivity on effective labor supply
Aprod=1; %TFP
RS = 0.9; %0.85; %important for fraction of entrepreneurs
theta = 0.33;
al = theta*RS; %capital share in production
be = (1-theta)*RS; %labor share in production
rstar=0.03; %fixed rental rate (small open economy)
discrate_workers = rstar; %Discount rate for workers
% New parameters for log(z) process - low persistence
Corr = 0.85;
nu = -log(Corr); %3;
%Var = 0.56^2/2;
Var = 0.3^2/(2*nu);
sig = 0.3;
sig2 = sig^2;
%sig2 = Var*2*nu; %0.041/0.0512*nu;
logzmean = 0;

%Grid points
J=30; %30; %Number of grid points z-dimension
I=200; %300; %Number of grid points a-dimension
T = 150; %150; %300; %Number of time periods
N = 150; %100; %300; %Number of steps in time dimension

% INITIAL DISTRIBUTION
contract=0.1; %Contraction of mean wealth by this amount
factor=0.1; % Scaling down the initial wealth 

% Vectors for flat tax
tau0_vec_tr = linspace(-0.4,-0.3,10)';
taubar_vec_tr = linspace(-0.4,-0.3,10)';
gal_vec_tr = [0];

% Vectors of parameters for optimal tax schedule
tax_vec_ss = linspace(0,0.1,20)'; %Testing for optimal steady state tax

tau0_vec = linspace(0.9,1.0,10)'; %Initial tax rate
taubar_vec = linspace(-0.8,-0.7,10)'; %Final tax rate
halflives_vec = [7 8 9];%[5 7 10]; %Tax schedule half life
% halflives_vec = [8]; %Tax schedule half life
% tau0_vec = linspace(-0.1,0,2)'; %Initial tax rate
% taubar_vec = linspace(0.1,0.2,2)'; %Final tax rate

gal_vec = -log(0.5)./halflives_vec;

%Time grid (uniform)
%timevec = linspace(0,T,N)';
%Time grid (non-uniform)
x = linspace(0,1,N)';
coeff = 50; power = 5; %Higher values give more unequal spacing of grid
xx  = x + coeff*x.^power;
xmax = max(xx); xmin = min(xx);
timevec = T/(xmax - xmin)*xx; %Non-uniform time vector
Deltavec = timevec(2:end)-timevec(1:end-1);
Deltavec = [Deltavec(1); Deltavec; Deltavec(end)];

%ITERATION PARAMETERS
flag=0; %Used to exit after non-convergence at any stage
maxit_price = 50; 
maxit_price_seq = 40;
maxit_vf=200;
maxit_dist=200;
crit_price = 10^(-6); %Convergence criterion steady state w
crit_price_seq = 10^(-5); %Convergence criterion sequence w
crit_dist = 10^(-8);
crit_vf = 10^(-8);
crit_transfer = 10^(-5);
crit_transfer_seq = 10^(-5);
Delta = 1000; %Time step (change for sequence problem)
%Search range for wage for initial and terminal problem
w0=1;
wmin0=1;
wmax0=5;

%Productivity grid (uniform)
zmin = 0.3;
zmax = 2.2;
z = linspace(zmin,zmax,J); %1xJ matrix
zz = ones(I,1)*z; %IxJ matrix of z repeated in I identical rows
dz = (zmax-zmin)/(J-1);
dzz = dz*ones(1,J);
dz2 = dz^2;
zs = reshape(zz,I*J,1);

%Asset grid (non-uniform)
eeps = 0.0001;
amin = eeps; %borrowing constraint
amax = 350; %130;
% a = linspace(amin,amax,I)'; %Ix1 matrix
% da = (amax-amin)/(I-1);
x = linspace(0,1,I)';
coeff = 5; power = 5; %15; %Higher values give more unequal spacing of grid
xx  = x + coeff*x.^power;
xmax = max(xx); xmin = min(xx);
a = (amax-amin)/(xmax - xmin)*xx + amin; %Non-uniform asset vector
aa = a*ones(1,J); %IxJ matrix of a repeated in J identical columns
daf = ones(I,1);
dab = ones(I,1);
daf(1:I-1) = a(2:I)-a(1:I-1);
dab(2:I) = a(2:I)-a(1:I-1);
daf(I)=daf(I-1); dab(1)=dab(2);
daaf = daf*ones(1,J); %Forward difference in grid
daab = dab*ones(1,J); %Backward difference
da_tilde = 0.5*(dab + daf); %Used for correct numerical integration
da_tilde(1) = 0.5*daf(1); da_tilde(I) = 0.5*dab(I);
da_stacked = reshape(da_tilde*ones(1,J),I*J,1);
grid_diag = spdiags(da_stacked,0,I*J,I*J);
as = reshape(aa,I*J,1);


%ORNSTEIN-UHLENBECK PROCESS dlog(zs) = -nu*log(z)dt + sig2*dW
%STATIONARY DISTRIBUTION IS log(z) ~ N(0,Var) WHERE Var = sig2/(2*nu)
%zmean = exp(logzmean + Var/2); %MEAN OF LOG-NORMAL DISTRIBUTION N(0,Var)
% Old parameters for log(z) process - high persistence
% nu = 0.02;
% sig2 = 0.041/0.0512*nu;
% logzmean = 0;

% Parameters for z process from log(z) process
mu = (nu*(logzmean - log(z)) + sig2/2).*z; %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*z.^2; %VARIANCE (FROM ITO'S LEMMA)


%CONSTRUCT MATRIX Bswitch SUMMARIZING EVOLUTION OF z
%yy = - s2/dz2 - mu/dz;
%chi =  s2/(2*dz2);    
%zeta = mu/dz + s2/(2*dz2);
chi =  - min(mu,0)/dz + s2/(2*dz2);
yy =  min(mu,0)/dz - max(mu,0)/dz - s2/dz2;
zeta = max(mu,0)/dz + s2/(2*dz2);

%This will be the upperdiagonal of the B_switch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)];
end

%This will be the center diagonal of the B_switch
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

%This will be the lower diagonal of the B_switch
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)];
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Bswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

%INITIALIZE MATRICES
Vaf = zeros(I,J);
Vab = zeros(I,J);
Vzf = zeros(I,J);
Vzb = zeros(I,J);
Vzz = zeros(I,J);

% PARETO DISTRIBUTED POISSON SHOCKS
if poissonshocks == 1
    % Distribution of shocks
    paretoparam = 1.1;
    shockcdf = @(x) (1-(zmin./x).^paretoparam)./(1-(zmin/zmax)^paretoparam); %Pareto cdf
    
%     zmin2 = zmin; %1/zmax;
%     zmax2 = 1.5;
%     shockcdf = @(x) myfunc(x,zmin2,zmax2);
    
    evalvec = zmin+([0 (linspace(1,J-1,J-1)-0.5) J-1])*dz;
    Ptil = shockcdf(evalvec(2:end))-shockcdf(evalvec(1:end-1));
    % Matrix for HJB algorithm
    aux=[Ptil; zeros(I-1,J)];
    aux2 = [zeros(1,I-1) aux(:)'];
    Cswitch = [];
    for i=1:I
        Cswitch = [Cswitch; aux2(I-i+1:I*J+I-i)];
    end
    Cswitch = poissonarrivalrate*repmat(Cswitch,J,1)-poissonarrivalrate*speye(I*J); %Just changed the signs here from v7 (since A changed sign)
    Bswitch = sparse(Cswitch+Bswitch);
end
CSwitch=Cswitch';