function outputs_noTVB = DAPM_noTVB(Rx,X1,X2,X3,X4,NOHET)

outputs_noTVB = [];

horsel = [1:120];
maxh = max(horsel);

%% Preliminary Calculations
Re = Rx(2:end,:);
[T,N] = size(Rx);
k1 = size(X1,2);
k2 = size(X2,2);
k3 = size(X3,2);
k4 = size(X4,2);
k = k1+k2+k3;
kP = k1+k2;
kF = k2+k3;
%%
X_full = [X1 X2 X3];
%%
X = X_full(2:end,:);
F = X(:,k1+1:end);
X_ = X_full(1:end-1,:);
mu_X_ = mean(X_)';
%%
Z_ = [ones(T-1,1) X_];
F_ = X_(:,k1+1:end);
Yps_FF = ([ones(T-1,1) F_]'*[ones(T-1,1) F_])/(T-1);
Yps_FF_inv = (T-1)*(([ones(T-1,1) F_]'*[ones(T-1,1) F_])\eye(kF+1));
X1_ = X_(:,1:k1);
Yps_11 = (X1_'*X1_)/(T-1);
Yps_F1 = ([ones(T-1,1) F_]'*X1_)/(T-1);
%% Estimation of VAR
Psi = (Z_'*Z_)\(Z_'*X);
PZ = Z_/(Z_'*Z_)*Z_'; 
MZ = eye(T-1) - PZ;
% Calculate VAR innovations
V = MZ*X;
U = V(:,1:k1+k2);
SigU = U'*U/(T-1);
SigV = V'*V/(T-1);
SigVU = V'*U/(T-1);
% Estimate SUR Model

W = [ones(T-1,1) F_ U];
W_bc = [W X1(1:end-1,:)];
A = (W'*W)\(W'*Re);
A_bc = (W_bc'*W_bc)\(W_bc'*Re);
E = Re - W*A;
E_bc = Re - W_bc*A_bc;
SigE = (E'*E)/(T-1);
B = A(end-kP+1:end,:);
B_bc = A_bc(end-k1-kP+1:end-k1,:);
%% Construct Robust Variance Estimator
Vrob = zeros(N*(kP+kF+1));
Vrob_diag = zeros(N*(kP+kF+1));
if ~NOHET,
    for t = 1:T-1,
        Vrob = Vrob + kron(W(t,:)'*W(t,:),E(t,:)'*E(t,:))/(T-1);
        Vrob_diag = Vrob_diag + kron(W(t,:)'*W(t,:),diag(diag(E(t,:)'*E(t,:))))/(T-1);
    end;
end;
Vrob = kron(inv(W'*W/T),eye(N))*Vrob*kron(inv(W'*W/T),eye(N))'; 
Vrob_diag = kron(inv(W'*W/T),eye(N))*Vrob_diag*kron(inv(W'*W/T),eye(N))'; 
% VB_hom = DO THESE LATER
% VB_hom_diag = 
VB_het = Vrob(end-N*kP+1:end,end-N*kP+1:end);
VB_het_diag = Vrob_diag(end-N*kP+1:end,end-N*kP+1:end);

%% Construct Robust Variance Estimator for Bias-Corrected Specification
Vrob_bc = zeros(N*(kP+kF+k1+1));
Vrob_bc_diag = zeros(N*(kP+kF+k1+1));
if ~NOHET,
    for t = 1:T-1,
        Vrob_bc = Vrob_bc + kron(W_bc(t,:)'*W_bc(t,:),E_bc(t,:)'*E_bc(t,:))/(T-1);
        Vrob_bc_diag = Vrob_bc_diag + kron(W_bc(t,:)'*W_bc(t,:),diag(diag(E_bc(t,:)'*E_bc(t,:))))/(T-1);
    end;
end;
Vrob_bc = kron(inv(W_bc'*W_bc/T),eye(N))*Vrob_bc*kron(inv(W_bc'*W_bc/T),eye(N))'; 
Vrob_bc_diag = kron(inv(W_bc'*W_bc/T),eye(N))*Vrob_bc_diag*kron(inv(W_bc'*W_bc/T),eye(N))'; 
% Note: We only require the upper left block of this matrix
Vrob_bc = Vrob_bc(1:N*(1+kP+kF),1:N*(1+kP+kF));
Vrob_bc_diag = Vrob_bc_diag(1:N*(1+kP+kF),1:N*(1+kP+kF));
% Calculate "OLS" Estimator of Lambda
Lam = (B*B')\(B*A(1:kF+1,:)');
F_tilde = [ones(T-1,1) F_]; 
Lamt = (F_tilde*Lam');
Premt = (F_tilde*Lam'*B);
% B_4S = (Lamt'*Lamt)\(Lamt'*Re);
% Calculate Bias-Corrected "OLS" Estimator of Lambda
Lam_bc = (B_bc*B_bc')\(B_bc*A_bc(1:kF+1,:)');
% Calculate "MLE" Estimator of Lambda and B
LL = A'*W'*W*A;
[evec, eigval] = eig(LL);
[blah, ind_evec] = sort(diag(eigval), 'descend');
L = evec(:,ind_evec(1:kP));
B_mle_tmp = L;
D_mle_tmp = B_mle_tmp'*A';
B_mle = B_mle_tmp*D_mle_tmp(:,kF+2:kF+1+kP);
B_mle = B_mle';
D_mle = D_mle_tmp(:,kF+2:kF+1+kP)\D_mle_tmp;
Lam_mle = D_mle(:,1:kF+1);
Lamt_mle = (F_tilde*Lam_mle');
Prem_mle = ([ones(T-1,1) F_]*Lam_mle'*B_mle);
% Calculate Bias-Corrected "MLE" Estimator of Lambda and B
B_mle_bc = [];
Lam_mle_bc = [];
if ~isempty(X1),
    M1 = T*(W'*W/T) - ((W'*X1(1:end-1,:)/T)/(X1(1:end-1,:)'*X1(1:end-1,:)/T))*(X1(1:end-1,:)'*W/T);
    LL_bc = A_bc(1:kP+kF+1,:)'*M1*A_bc(1:kP+kF+1,:);
    [evec_bc, eval_bc] = eig(LL_bc);
    [blah, ind_evec_bc] = sort(diag(eval_bc), 'descend');
    L_bc = evec_bc(:,ind_evec_bc(1:kP));
    B_mle_bc_tmp = L_bc;
    D_mle_bc_tmp = B_mle_bc_tmp'*A_bc';
    B_mle_bc = B_mle_bc_tmp*D_mle_bc_tmp(:,kF+2:kF+1+kP);
    B_mle_bc = B_mle_bc';
    D_mle_bc = D_mle_bc_tmp(:,kF+2:kF+1+kP)\D_mle_bc_tmp;
    Lam_mle_bc = D_mle_bc(:,1:kF+1);
    Lamt_mle_bc = (F_tilde*Lam_mle_bc');
end

%% Estimate with X4 Variable
B_4 = [];
Lam_4 = [];
B_mle_4 = [];
Lam_mle_4=[];
if ~isempty(X4)
    W_4 = [ones(T-1,1) F_ U X4(2:end,:)];
    A_4 = (W_4'*W_4)\(W_4'*Re);
    B_4 = A_4(end-k4-kP+1:end-k4,:);
    Lam_4 = (B_4*B_4')\(B_4*A_4(1:kF+1,:)');
    %%
    M4 = T*(W'*W/T) - ((W'*X4(2:end,:)/T)/(X4(2:end,:)'*X4(2:end,:)/T))*(X4(2:end,:)'*W/T);
    LL_4 = A_4(1:kP+kF+1,:)'*M4*A_4(1:kP+kF+1,:);
    [evec_4, eval_4] = eig(LL_4);
    [blah, ind_evec_4] = sort(diag(eval_4), 'descend');
    L_4 = evec_4(:,ind_evec_4(1:kP));
    B_mle_4_tmp = L_4;
    D_mle_4_tmp = B_mle_4_tmp'*A_4';
    B_mle_4 = B_mle_4_tmp*D_mle_4_tmp(:,kF+2:kF+1+kP);
    B_mle_4 = B_mle_4';
    D_mle_4 = D_mle_4_tmp(:,kF+2:kF+1+kP)\D_mle_4_tmp;
    Lam_mle_4 = D_mle_4(:,1:kF+1);
end;

%% Calculate Standard Errors/t-Statistics for "OLS" Estimator of Lambda
VLam_hom = kron(Yps_FF_inv,SigU) + kron((Yps_FF_inv+Lam'/SigU*Lam),(B*B')\B*SigE*B'/(B*B'));
VLam_hom_diag = kron(Yps_FF_inv,SigU) + kron((Yps_FF_inv+Lam'/SigU*Lam),(B*B')\B*diag(diag(SigE))*B'/(B*B'));
seLam_hom = (reshape(diag(VLam_hom./(T-1)),k1+k2,k2+k3+1).^(1/2));
seLam_hom_diag = (reshape(diag(VLam_hom_diag./(T-1)),k1+k2,k2+k3+1).^(1/2));
tLam_hom = Lam./seLam_hom;
tLam_hom_diag = Lam./seLam_hom_diag;
Lamt_hom_upper = []; Lamt_hom_lower = [];
if kF>0,
    for t = 1:T-1,
        seLamt_hom(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_hom*kron(F_tilde(t,:),eye(kP))'./(T-1)));    
        seLamt_hom_diag(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_hom_diag*kron(F_tilde(t,:),eye(kP))'./(T-1)));
    end;
    Lamt_hom_upper = Lamt+1.96*seLamt_hom; Lamt_hom_lower = Lamt-1.96*seLamt_hom;
    Lamt_hom_diag_upper = Lamt+1.96*seLamt_hom_diag; Lamt_hom_lower = Lamt-1.96*seLamt_hom_diag;
end;
%% Calculate Standard Errors/t-Statistics for "OLS" Estimator of B
seB_het = reshape(diag(VB_het)./(T-1),kP,N).^(1/2);
seB_het_diag = reshape(diag(VB_het_diag)./(T-1),kP,N).^(1/2);
tB_het = B./seB_het;
tB_het_diag = B./seB_het_diag;
%% confidence bands for lambdat
H = [kron(eye(kF+1),(B*B')\B) -kron(Lam',(B*B')\B)];
VLam_het = kron(Yps_FF_inv,SigU) + H*Vrob*H';
VLam_het_diag = kron(Yps_FF_inv,SigU) + H*Vrob_diag*H';
Lamt_het_upper = []; Lamt_het_lower = [];
Lamt_het_diag_upper = []; Lamt_het_diag_lower = [];
if kF > 0
    for t = 1:T-1,
        seLamt_het(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_het*kron(F_tilde(t,:),eye(kP))'./(T-1)));
        seLamt_het_diag(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_het_diag*kron(F_tilde(t,:),eye(kP))'./(T-1)));
    end;
    Lamt_het_upper = Lamt+1.96*seLamt_het; Lamt_het_lower = Lamt-1.96*seLamt_het;
    Lamt_het_diag_upper = Lamt+1.96*seLamt_het_diag; Lamt_het_diag_lower = Lamt-1.96*seLamt_het_diag;
end;
seLam_het = (reshape(diag(VLam_het./(T-1)),k1+k2,k2+k3+1).^(1/2));
seLam_het_diag = (reshape(diag(VLam_het_diag./(T-1)),k1+k2,k2+k3+1).^(1/2));
tLam_het = Lam./seLam_het;
tLam_het_diag = Lam./seLam_het_diag;


%% Calculate Standard Errors/t-Statistics for Bias-Corrected "OLS" Estimator of Lambda
H_bc = [kron(eye(kF+1),(B_bc*B_bc')\B_bc) -kron(Lam_bc',(B_bc*B_bc')\B_bc)];
VLam_bc_het = kron(eye(kF+1)/(Yps_FF - Yps_F1/Yps_11*Yps_F1'),SigU) + H_bc*Vrob_bc*H_bc';
VLam_bc_het_diag = kron(eye(kF+1)/(Yps_FF - Yps_F1/Yps_11*Yps_F1'),SigU) + H_bc*Vrob_bc_diag*H_bc';
Lamt_bc_het_upper = []; Lamt_bc_het_lower = [];
Lamt_bc_het_diag_upper = []; Lamt_bc_het_diag_lower = [];
if kF > 0,
    for t = 1:T-1,
        seLamt_bc_het(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_bc_het*kron(F_tilde(t,:),eye(kP))'./(T-1)));
        seLamt_bc_het_diag(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_bc_het_diag*kron(F_tilde(t,:),eye(kP))'./(T-1)));
    end;
    Lamt_bc_het_upper = Lamt+1.96*seLamt_bc_het; Lamt_bc_het_lower = Lamt-1.96*seLamt_bc_het;
    Lamt_bc_het_diag_upper = Lamt+1.96*seLamt_bc_het_diag; Lamt_bc_het_diag_lower = Lamt-1.96*seLamt_bc_het_diag;
end;
seLam_bc_het = (reshape(diag(VLam_bc_het./(T-1)),k1+k2,k2+k3+1).^(1/2));
seLam_bc_het_diag = (reshape(diag(VLam_bc_het_diag./(T-1)),k1+k2,k2+k3+1).^(1/2));
tLam_bc_het = Lam_bc./seLam_bc_het;
tLam_bc_het_diag = Lam_bc./seLam_bc_het_diag;
%% Calculate Standard Errors/t-Statistics for "MLE" Estimator of Lambda
VLam_mle_hom = kron(Yps_FF_inv,SigU) + kron((Yps_FF_inv+Lam_mle'/SigU*Lam_mle),(B_mle*B_mle')\B_mle*SigE*B_mle'/(B_mle*B_mle'));
VLam_mle_hom_diag = kron(Yps_FF_inv,SigU) + kron((Yps_FF_inv+Lam_mle'/SigU*Lam_mle),(B_mle*B_mle')\B_mle*diag(diag(SigE))*B_mle'/(B_mle*B_mle'));
Lamt_mle_hom_upper = []; Lamt_mle_hom_lower = [];
Lamt_mle_hom_diag_upper = []; Lamt_mle_hom_diag_lower = [];
if kF > 0
    for t = 1:T-1,
        seLamt_mle_hom(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_mle_hom*kron(F_tilde(t,:),eye(kP))'./(T-1)));
        seLamt_mle_hom_diag(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_mle_hom_diag*kron(F_tilde(t,:),eye(kP))'./(T-1)));
    end;
    Lamt_mle_hom_upper = Lamt_mle+1.96*seLamt_mle_hom; Lamt_mle_hom_lower = Lamt_mle-1.96*seLamt_mle_hom;
    Lamt_mle_hom_diag_upper = Lamt_mle+1.96*seLamt_mle_hom_diag; Lamt_mle_hom_diag_lower = Lamt_mle-1.96*seLamt_mle_hom_diag;
end;
seLam_mle_hom = (reshape(diag(VLam_mle_hom./(T-1)),k1+k2,k2+k3+1).^(1/2));
seLam_mle_hom_diag = (reshape(diag(VLam_mle_hom_diag./(T-1)),k1+k2,k2+k3+1).^(1/2));
tLam_mle_hom = Lam_mle./seLam_mle_hom;
tLam_mle_hom_diag = Lam_mle./seLam_mle_hom_diag;
%%
H_mle = [kron(eye(kF+1),(B_mle*B_mle')\B_mle) -kron(Lam_mle',(B_mle*B_mle')\B_mle)];
VLam_mle_het = kron(Yps_FF_inv,SigU) + H_mle*Vrob*H_mle';
VLam_mle_het_diag = kron(Yps_FF_inv,SigU) + H_mle*Vrob_diag*H_mle';
Lamt_mle_het_upper = []; Lamt_mle_het_lower = [];
Lamt_mle_het_diag_upper = []; Lamt_mle_het_diag_lower = [];
if kF > 0,
    for t = 1:T-1,
        seLamt_mle_het(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_mle_het*kron(F_tilde(t,:),eye(kP))'./(T-1)));
        seLamt_mle_het_diag(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_mle_het_diag*kron(F_tilde(t,:),eye(kP))'./(T-1)));
    end;
    Lamt_mle_het_upper = Lamt_mle+1.96*seLamt_mle_het; Lamt_mle_het_lower = Lamt_mle-1.96*seLamt_mle_het;
    Lamt_mle_het_diag_upper = Lamt_mle+1.96*seLamt_mle_het_diag; Lamt_mle_het_diag_lower = Lamt_mle-1.96*seLamt_mle_het_diag;
end;
seLam_mle_het = (reshape(diag(VLam_mle_het./(T-1)),k1+k2,k2+k3+1).^(1/2));
seLam_mle_het_diag = (reshape(diag(VLam_mle_het_diag./(T-1)),k1+k2,k2+k3+1).^(1/2));
tLam_mle_het = Lam_mle./seLam_mle_het;
tLam_mle_het_diag = Lam_mle./seLam_mle_het_diag;
%% Variance of Beta for MLE estimator
LUL_inv = eye(kP)/(Lam_mle*Yps_FF*Lam_mle' + SigU);
HB_mle = kron([LUL_inv*Lam_mle LUL_inv]*(W'*W/T),eye(N)) - kron(LUL_inv*Lam_mle*Yps_FF, B_mle')*H_mle;
VB_mle_het = HB_mle*Vrob*HB_mle';
VB_mle_het_diag = HB_mle*Vrob_diag*HB_mle';
%% Calculate Standard Errors/t-Statistics for "MLE" Estimator of B
seB_mle_het = reshape(diag(VB_mle_het)./(T-1),kP,N).^(1/2);
seB_mle_het_diag = reshape(diag(VB_mle_het_diag)./(T-1),kP,N).^(1/2);
tB_mle_het = B_mle./seB_mle_het;
tB_mle_het_diag = B_mle./seB_mle_het_diag;
%% % Wald test for an entire row of B being jointly zero
k12 = k1 + k2;
for j = 1:k12,
    B_vec = B(j,:)';
    B_indx = [j:k12:N*k12];                
    WB_het(j) = (T-1)*B_vec'/VB_het(B_indx,B_indx)*B_vec;
    pWB_het(j) = 1-chi2cdf(WB_het(j),length(B_vec));                
    WB_het_diag(j) = (T-1)*B_vec'/VB_het_diag(B_indx,B_indx)*B_vec;
    pWB_het_diag(j) = 1-chi2cdf(WB_het_diag(j),length(B_vec));                

    B_mle_vec = B_mle(j,:)';
    WB_mle_het(j) = (T-1)*B_mle_vec'/VB_mle_het(B_indx,B_indx)*B_mle_vec;
    pWB_mle_het(j) = 1-chi2cdf(WB_mle_het(j),length(B_mle_vec));                
    WB_mle_het_diag(j) = (T-1)*B_mle_vec'/VB_mle_het_diag(B_indx,B_indx)*B_mle_vec;
    pWB_mle_het_diag(j) = 1-chi2cdf(WB_mle_het_diag(j),length(B_mle_vec));                
end;

%% Calculate Standard Errors/t-Statistics for Bias-Corrected "MLE" Estimator of Lambda
if ~isempty(X1),
    H_mle_bc = [kron(eye(kF+1),(B_mle_bc*B_mle_bc')\B_mle_bc) -kron(Lam_mle_bc',(B_mle_bc*B_mle_bc')\B_mle_bc)];
    VLam_mle_bc_het = kron(eye(kF+1)/(Yps_FF - Yps_F1/Yps_11*Yps_F1'),SigU) + H_mle_bc*Vrob_bc*H_mle_bc';
    VLam_mle_bc_het_diag = kron(eye(kF+1)/(Yps_FF - Yps_F1/Yps_11*Yps_F1'),SigU) + H_mle_bc*Vrob_bc_diag*H_mle_bc';
    Lamt_mle_bc_het_upper = []; Lamt_mle_bc_het_lower = [];
    Lamt_mle_bc_het_diag_upper = []; Lamt_mle_bc_het_diag_lower = [];
    if kF>0,
        for t = 1:T-1,
            seLamt_mle_bc_het(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_mle_bc_het*kron(F_tilde(t,:),eye(kP))'./(T-1)));
            seLamt_mle_bc_het_diag(t,:) = sqrt(diag(kron(F_tilde(t,:),eye(kP))*VLam_mle_bc_het_diag*kron(F_tilde(t,:),eye(kP))'./(T-1)));
        end;
    end;
    seLam_mle_bc_het = (reshape(diag(VLam_mle_bc_het./(T-1)),k1+k2,k2+k3+1).^(1/2));
    seLam_mle_bc_het_diag = (reshape(diag(VLam_mle_bc_het_diag./(T-1)),k1+k2,k2+k3+1).^(1/2));
    tLam_mle_bc_het = Lam_mle_bc./seLam_mle_bc_het;
    tLam_mle_bc_het_diag = Lam_mle_bc./seLam_mle_bc_het_diag;

    % Bias-Corrected Versions
    Lam1_bc_tilde = [zeros(kP,k1) Lam_bc(:,2:end)];
    mu_F_tilde = [1 mu_X_(k1+1:end)']';
    mu_1 = mu_X_(1:k1);
    CovLam_bar_bc_hat1 = (mu_F_tilde'/(Yps_FF - Yps_F1/Yps_11*Yps_F1')*mu_F_tilde)*Lam1_bc_tilde/(eye(k)-Psi(2:end,:)')*SigVU;
    CovLam_bar_bc_hat2 = -(mu_1'/Yps_11*Yps_F1'/(Yps_FF - Yps_F1/Yps_11*Yps_F1')*mu_F_tilde)*Lam1_bc_tilde/(eye(k)-Psi(2:end,:)')*SigVU;
    VLam_bar_bc_het = kron(mu_F_tilde',eye(kP))*VLam_bc_het*kron(mu_F_tilde',eye(kP))' + Lam1_bc_tilde/(eye(k)-Psi(2:end,:)')*SigV/((eye(k)-Psi(2:end,:)')')*Lam1_bc_tilde' + ...
        CovLam_bar_bc_hat1 + CovLam_bar_bc_hat1' + CovLam_bar_bc_hat2 + CovLam_bar_bc_hat2';                
    seLam_bar_bc_het = ((diag(VLam_bar_bc_het)./(T-1)).^(1/2));
    Lam_bar_bc = (Lam_bc(:,1) + Lam1_bc_tilde*mu_X_);
    tLam_bar_bc_het = Lam_bar_bc./seLam_bar_bc_het;

    Lam1_mle_bc_tilde = [zeros(kP,k1) Lam_mle_bc(:,2:end)];
    mu_F_tilde = [1 mu_X_(k1+1:end)']';
    mu_1 = mu_X_(1:k1);
    CovLam_bar_mle_bc_hat1 = (mu_F_tilde'/(Yps_FF - Yps_F1/Yps_11*Yps_F1')*mu_F_tilde)*Lam1_mle_bc_tilde/(eye(k)-Psi(2:end,:)')*SigVU;
    CovLam_bar_mle_bc_hat2 = -(mu_1'/Yps_11*Yps_F1'/(Yps_FF - Yps_F1/Yps_11*Yps_F1')*mu_F_tilde)*Lam1_mle_bc_tilde/(eye(k)-Psi(2:end,:)')*SigVU;
    VLam_bar_mle_bc_het = kron(mu_F_tilde',eye(kP))*VLam_mle_bc_het*kron(mu_F_tilde',eye(kP))' + Lam1_mle_bc_tilde/(eye(k)-Psi(2:end,:)')*SigV/((eye(k)-Psi(2:end,:)')')*Lam1_mle_bc_tilde' + ...
        CovLam_bar_mle_bc_hat1 + CovLam_bar_mle_bc_hat1' + CovLam_bar_mle_bc_hat2 + CovLam_bar_mle_bc_hat2';                
    seLam_bar_mle_bc_het = ((diag(VLam_bar_mle_bc_het)./(T-1)).^(1/2));
    Lam_bar_mle_bc = (Lam_mle_bc(:,1) + Lam1_mle_bc_tilde*mu_X_);
    tLam_bar_mle_bc_het = Lam_bar_mle_bc./seLam_bar_mle_bc_het;

end;

%% Wald test for rows of Lambda being jointly zero
vartype = {'het','het_diag','hom','hom_diag'};
esttype = {'Lam','Lam_bc','Lam_mle','Lam_mle_bc'};
for i = 1:length(esttype),
    for j = 1:length(vartype),
        if ~isempty(strfind(esttype{i},'bc')) && ~isempty(strfind(vartype{j},'hom')),
            continue;
        end;
        for l = 1:kP,            
            dum = zeros(size(Lam));
            dum(l,:) = ones(1,size(Lam,2));
            dum1 = zeros(size(Lam));
            dum1(l,2:end) = ones(1,size(Lam,2)-1);
            lindx = find(dum);
            l1indx = find(dum1);
            eval(['Wald_',esttype{i},'_',vartype{j},'{l} = (T-1)*Lam(l,:)/V',esttype{i},'_',vartype{j},'(lindx,lindx)*Lam(l,:)',''';']);    
            eval(['pWald_',esttype{i},'_',vartype{j},'{l} = 1 - chi2cdf(Wald_',esttype{i},'_',vartype{j},'{l},size(Lam,2));']);    
            eval(['Wald_',strrep(esttype{i},'Lam','Lam1'),'_',vartype{j},'{l} = (T-1)*Lam(l,2:end)/V',esttype{i},'_',vartype{j},'(l1indx,l1indx)*Lam(l,2:end)',''';']);    
            eval(['pWald_',strrep(esttype{i},'Lam','Lam1'),'_',vartype{j},'{l} = 1 - chi2cdf(Wald_',strrep(esttype{i},'Lam','Lam1'),'_',vartype{j},'{l},size(Lam,2)-1);']);    
%             for c = 2:kF+1,
%                 dum = zeros(size(Lam));
%                 dum(l,c) = 1;            
%                 lcindx = find(dum);
%                 eval(['Wald_',strrep(esttype{i},'Lam','Lam11'),'_',vartype{j},'{l} = (T-1)*Lam(l,c)/V',esttype{i},'_',vartype{j},'(lcindx,lcindx)*Lam(l,c)','''']);    
%                 eval(['pWald_',strrep(esttype{i},'Lam','Lam11'),'_',vartype{j},'{l} = 1 - chi2cdf(Wald_',strrep(esttype{i},'Lam','Lam11'),'_',vartype{j},'{l},1)']);    
%             end;
        end;
    end;
end;

%% Compute test statistic for Lambda_bar
Lam1_tilde = [zeros(kP,k1) Lam(:,2:end)];
mu_F_tilde = [1 mu_X_(k1+1:end)']';
CovLam_bar_hat = (mu_F_tilde'*Yps_FF_inv*mu_F_tilde)*Lam1_tilde/(eye(k)-Psi(2:end,:)')*SigVU;
VLam_bar_het = kron(mu_F_tilde',eye(kP))*VLam_het*kron(mu_F_tilde',eye(kP))' + Lam1_tilde/(eye(k)-Psi(2:end,:)')*SigV/((eye(k)-Psi(2:end,:)')')*Lam1_tilde' + CovLam_bar_hat + CovLam_bar_hat';                
seLam_bar_het = ((diag(VLam_bar_het)./(T-1)).^(1/2));
Lam_bar = (Lam(:,1) + Lam1_tilde*mu_X_);
tLam_bar_het = Lam_bar./seLam_bar_het;

Lam1_mle_tilde = [zeros(kP,k1) Lam_mle(:,2:end)];
mu_F_tilde = [1 mu_X_(k1+1:end)']';
CovLam_bar_mle_hat = (mu_F_tilde'*Yps_FF_inv*mu_F_tilde)*Lam1_mle_tilde/(eye(k)-Psi(2:end,:)')*SigVU;
VLam_bar_mle_het = kron(mu_F_tilde',eye(kP))*VLam_mle_het*kron(mu_F_tilde',eye(kP))' + Lam1_mle_tilde/(eye(k)-Psi(2:end,:)')*SigV/((eye(k)-Psi(2:end,:)')')*Lam1_mle_tilde' + CovLam_bar_mle_hat + CovLam_bar_mle_hat';                
seLam_bar_mle_het = ((diag(VLam_bar_mle_het)./(T-1)).^(1/2));
Lam_bar_mle = (Lam_mle(:,1) + Lam1_mle_tilde*mu_X_);
tLam_bar_mle_het = Lam_bar_mle./seLam_bar_mle_het;

%% Calculate Fitted Values
% Rxhat = ([ones(T-1,1) F_]*Lam'+ U)*B;
Rxhat_mle = ([ones(T-1,1) F_]*Lam_mle'+ U)*B_mle;
if isempty([X1 X3]), 
    Rxhat = ([ones(T-1,1) F_]*Psi+ U)*B;
else
    Rxhat = ([ones(T-1,1) F_]*Lam'+ U)*B;
end;

Rxhat_bc = [];
Rxhat_mle_bc = [];
E2_bc = [];
E2_mle_bc = [];
if ~isempty(X1),
    Rxhat_bc = ([ones(T-1,1) F_]*Lam_bc'+ U)*B_bc + X1(1:end-1,:)*(A_bc(end-k1+1:end,:) + (X1(1:end-1,:)'*X1(1:end-1,:)\(X1(1:end-1,:)'*W)*(A_bc(1:end-k1,:)-[Lam_bc eye(kP)]'*B_bc)));
    Rxhat_mle_bc = ([ones(T-1,1) F_]*Lam_mle_bc'+ U)*B_mle_bc + X1(1:end-1,:)*(A_bc(end-k1+1:end,:) + (X1(1:end-1,:)'*X1(1:end-1,:)\(X1(1:end-1,:)'*W)*(A_bc(1:end-k1,:)-[Lam_mle_bc eye(kP)]'*B_mle_bc)));
    %Rxhat_bc = ([ones(T-1,1) F_]*Lam_bc'+ U)*B_bc + X1(1:end-1,:)*A_bc(end-k1+1:end,:);
    E2_bc = Re - Rxhat_bc;
    E2_mle_bc = Re - Rxhat_mle_bc;
end;
E2 = Re - Rxhat;
E2_mle = Re - Rxhat_mle;
%
Rxhat_4 = [];
Rxhat_mle_4 = [];
if ~isempty(X4)
    Rxhat_4 = ([ones(T-1,1) F_]*Lam_4'+ U)*B_4 + X4(2:end,:)*(A_4(end-k4+1:end,:) + (X4(2:end,:)'*X4(2:end,:)\(X4(2:end,:)'*W)*(A_4(1:end-k4,:)-[Lam_4 eye(kP)]'*B_4)));    
    Rxhat_mle_4 = ([ones(T-1,1) F_]*Lam_mle_4'+ U)*B_mle_4 + X4(2:end,:)*(A_4(end-k4+1:end,:) + (X4(2:end,:)'*X4(2:end,:)\(X4(2:end,:)'*W)*(A_4(1:end-k4,:)-[Lam_mle_4 eye(kP)]'*B_mle_4)));    
end

%% Principal components of returns and fitted return errors
[LoadingsRx,PCsRx,evalsRx] = princomp(demean(Rx));
sharePCRx = evalsRx(1:5)./sum(evalsRx(1:5));
PCsRxTotalVar = sum(evalsRx);

[LoadingsE_mle,PCsE_mle,evalsE_mle] = princomp(demean(Re-Rxhat_mle));
sharePCE_mle = evalsE_mle(1:5)./sum(evalsE_mle(1:5));


[LoadingsE,PCsE,evalsE] = princomp(demean(E));
sharePCE = evalsE(1:5)./sum(evalsE(1:5));
PCsETotalVar = sum(evalsE);

LR = T*(sum(log(evalsRx)-log(evalsE)));
pLR = 1-chi2cdf(LR,3*N);  
LRcv = chi2inv(.95, [1*N 2*N 3*N]);

% %% Fama-MacBeth estimator and variance
% V_FM = X_full - repmat(mean(X),T,1);
% W_FM = [ones(T,1) V_FM];
% A_FM = (W_FM'*W_FM)\(W_FM'*Rx);
% E_FM = Rx - W_FM*A_FM;
% SigE_FM = (E_FM'*E_FM)/T;
% B_FM = A_FM(end-kP+1:end,:);
% SigV = V_FM'*V_FM/T;
% Vrob_FM_hom = kron(inv(W_FM'*W_FM/T),SigE_FM);
% % Construct Robust Variance Estimator for FM estimator
% Vrob_FM = zeros(N*size(W_FM,2));
% Vrob_FM_diag = zeros(N*size(W_FM,2));
% if ~NOHET,
%     for t = 1:T-1,
%         Vrob_FM = Vrob_FM + kron(W_FM(t,:)'*W_FM(t,:),E_FM(t,:)'*E_FM(t,:))/(T-1);
%         Vrob_FM_diag = Vrob_FM_diag + kron(W_FM(t,:)'*W_FM(t,:),diag(diag(E_FM(t,:)'*E_FM(t,:))))/(T-1);
%     end;
% end;
% Vrob_FM_het = kron(inv(W_FM'*W_FM/T),eye(N))*Vrob_FM*kron(inv(W_FM'*W_FM/T),eye(N))'; 
% Lam_FM = (B_FM*B_FM')\(B_FM*A_FM(1,:)');
% %Lam_FM = (B_FM*B_FM')\(B_FM*mean(Re)');
% H_FM = [(B_FM*B_FM')\B_FM -kron(Lam_FM',(B_FM*B_FM')\B_FM)];
% VLam_FM_het = SigV + H_FM*Vrob_FM_het*H_FM';
% seLam_FM_het = (reshape(diag(VLam_FM_het./(T-1)),k1+k2,1).^(1/2));
% tLam_FM_het = Lam_FM./seLam_FM_het;
% 
% VLam_FM_hom = SigV + H_FM*Vrob_FM_hom*H_FM';
% VLam_FM_hom2 = SigV + kron(Lam_FM'*inv(SigV)*Lam_FM + 1,(B_FM*B_FM')\B_FM*SigE_FM*B_FM'*(B_FM*B_FM')); 
% seLam_FM_hom = (reshape(diag(VLam_FM_hom./(T-1)),k1+k2,1).^(1/2));
% tLam_FM_hom = Lam_FM./seLam_FM_hom;

%% One-Step Ahead Forecasts
F_end = [];
if ~isempty([X1 X2]) 
    F_end = F(end,:);
end
% Rxfcst = [1 F_end]*Lam'*B;
% Rxfcst_mle = [1 F_end]*Lam_mle'*B_mle;
% if ~isempty(X1)
%     Rxfcst_bc = [1 F_end]*Lam_bc'*B_bc + X1(end,:)*(A_bc(end-k1+1:end,:) + (X1(1:end-1,:)'*X1(1:end-1,:)\(X1(1:end-1,:)'*W)*(A_bc(1:end-k1,:)-[Lam_bc eye(kP)]'*B_bc)));
%     Rxfcst_mle_bc = [1 F_end]*Lam_mle_bc'*B_mle_bc + X1(end,:)*(A_bc(end-k1+1:end,:) + (X1(1:end-1,:)'*X1(1:end-1,:)\(X1(1:end-1,:)'*W)*(A_bc(1:end-k1,:)-[Lam_mle_bc eye(kP)]'*B_mle_bc)));
% end;
Rxfcst = [1 mean(F)]*Lam'*B;
Rxfcst_mle = [1 mean(F)]*Lam_mle'*B_mle;
if ~isempty(X1)
    Rxfcst_bc = [1 mean(F)]*Lam_bc'*B_bc + mean(X1)*(A_bc(end-k1+1:end,:) + (X1(1:end-1,:)'*X1(1:end-1,:)\(X1(1:end-1,:)'*W)*(A_bc(1:end-k1,:)-[Lam_bc eye(kP)]'*B_bc)));
    Rxfcst_mle_bc = [1 mean(F)]*Lam_mle_bc'*B_mle_bc + mean(X1)*(A_bc(end-k1+1:end,:) + (X1(1:end-1,:)'*X1(1:end-1,:)\(X1(1:end-1,:)'*W)*(A_bc(1:end-k1,:)-[Lam_mle_bc eye(kP)]'*B_mle_bc)));
end;

% Compute twelve-period ahead expected excess returns for all test assets
X_fcsts = zeros(T-1,k,maxh);
F_tilde_fcsts = zeros(T-1,1+kF,maxh);
mu = Psi(1,:)'; Phi = Psi(2:end,:)';
mu_fcsts = zeros(k,maxh);
mu_fcsts(:,1) = mu;
Phi_fcsts = zeros(k,k,maxh);
for h = 1:maxh,
    if h > 2,
        mu_fcsts(:,h) = [eye(k) + sum(Phi_fcsts(:,:,1:h-1),3)]*mu;
    end;
    Phi_fcsts(:,:,h) = Phi^h;
    X_fcsts(:,:,h) = Z_*[mu_fcsts(:,h) Phi_fcsts(:,:,h)]';
    F_tilde_fcsts(:,:,h) = [ones(T-1,1) X_fcsts(:,k1+1:end,h)];    
end;    

for hor = 1:length(horsel),
    h = horsel(hor);
    Premth{hor} = ones(T-1,N); Rxh{hor} = ones(T-1,N); 
    Premth_MKT{hor} = ones(T-1,1);
    for hh = 1:h,
        Lamt_h = F_tilde_fcsts(:,:,hh)*Lam';
        Premth_MKT{h} = Premth_MKT{h}.*(ones(T-1,1) + Lamt_h(:,1));
    end;    
    Premth_MKT{hor} = 100*(Premth_MKT{hor}.^(12/h)-1); % produce annual rates of return
    % These really are forward expected excess returns
    for i = 1:N,
        Premth{hor}(:,i) = Premth{hor}(:,i).*(ones(T-1,1) + sum(repmat(B(:,i),1,T-1)'.*(F_tilde_fcsts(:,:,h)*Lam'),2));        
        Rxh{hor}(:,i) = Rxh{hor}(:,i).*(ones(T-1,1) + [Re(h:end,i); NaN(h-1,1)]);        
    end; % for i
    Premth{hor} = 100*(Premth{hor}-1);
    Rxh{hor} = 100*(Rxh{hor}-1);
end; % for hor
    
    
%% Outputs
varlist = char(who);
for i = 1:size(varlist,1),
    eval([deblank(varlist(i,:)),'_CB = ',deblank(varlist(i,:)),';']);
    eval(['outputs_noTVB.',deblank(varlist(i,:)),'_CB = ',deblank(varlist(i,:)),'_CB;']);    
end;

