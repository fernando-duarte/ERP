% This example replicates the results of Ang and Kristensen (2011),
% "Testing Conditional Factor Models"
clear

%% load the data

% The variable X has the Fama-French factors. The
% variable Y has decile portfolios sorted on book-to-market. All series are
% at a daily frequency, July 1963 to December 2007. Data obtained in
% January 2012, which is slightly different from Ang and Kristensen (CRSP
% data changes over time)
load example.mat

%% Clean the data

Y = Y - repmat(rf,1,size(Y,2)); % make excess return
Y = Y/100; 
X = X/100;
keep = ~ ((year<=1964 & month<=6) | year == 1963 | year == 2007);
Ysub = Y(keep,:); % eliminate first and last year
Xsub = X(keep,:);
%% Reproduce Table 1
% Panel A
[mean(Xsub)'*252 std(Xsub)'*sqrt(252)]
corr(Xsub)

% Panel B
[mean(Ysub)'*252 std(Ysub)'*sqrt(252)]
for i = 1:10
    coeff(i,:) = regress(Ysub(:,i)*252,[ones(size(Xsub,1),1) Xsub(:,1)*252]);
end
coeff
%% Find optimal bandwidth
[hAlphaCAPM hBetaCAPM hLongRunAlphaCAPM hLongRunBetaCAPM] = optimal_bandwidth(Ysub(:,1),Xsub(:,1:2),'one-sided exponential');

[hAlphaCAPM hBetaCAPM hLongRunAlphaCAPM hLongRunBetaCAPM] = optimal_bandwidth(Ysub,Xsub(:,1),'one-sided-uniform');
[hAlphaCAPM1 hBetaCAPM1 hLongRunAlphaCAPM1 hLongRunBetaCAPM1] = optimal_bandwidth(Y,X(:,1));
[hAlphaCAPM2 hBetaCAPM2 hLongRunAlphaCAPM2 hLongRunBetaCAPM2] = optimal_bandwidth(Y(:,1)*252,X(:,1)*252,'one-sided-uniform',2);

[hAlphaFF hBetaFF hLongRunAlphaFF hLongRunBetaFF] = optimal_bandwidth(Y,X);
[hAlphaFF1 hBetaFF1 hLongRunAlphaFF1 hLongRunBetaFF1] = optimal_bandwidth(Ysub,Xsub);

 %% Find conditional coefficients
 bandwidth = [.0474 .0989 0.0349 0.0294 0.0379 0.0213 0.0188 0.0213 0.0160 0.0182];
 [alphaCAPM betaCAPM statsCAPM] = rolling_regress(Ysub,Xsub(:,1:2),'one-sided-uniform',hBetaCAPM2);
 [alphaFF betaFF statsFF] = rolling_regress(Y,X,'gaussian',bandwidth);
 
 %% Run test
test = 'Conditional Alpha';
[pvalue point_estimate conf_intvl st] = hypothesis_test(alphaCAPM,betaCAPM,statsCAPM,test,0.05);
test = 'Conditional Beta';
[pvalue point_estimate conf_intvl st] = hypothesis_test(alphaCAPM,betaCAPM,statsCAPM,test,0.05);

%% Plots
[T,J] = size(betaCAPM{1});
plot(1:T,betaCAPM{1},1:T,cell2mat(conf_intvl{1}))