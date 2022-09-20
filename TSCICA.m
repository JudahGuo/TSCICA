function [IC, w] = TSCICA(X, ref,reft,threshold, w0, mu1,mu2, lambda, maxIter, OverValue)
% Constrained ICA for extracting one source signal using the reference signal.
%
% Command:
%   [y, w] = cICA(whitesig, ref,reft,threshold, w0, mu1,mu2, lambda, maxIter, OverValue)
%
% Parameters:
%            IC --- extracted source signal, according to the given optimal time delay.
%            w --- corresponding weight vector, that is, y=w'*X.
%            X --- prewhitened observed mixed signals. Each row is a mixed signal.     
%            ref ---spatial reference signal
%            reft ---temporal reference signal
%    threshold --- the closeness of the desired signal and the reference signal is less than 
%                  or equal to the threshold, i.e., closeness(y,ref) - threshold <= 0
%           w0 --- the initial weight vector
% learningRate --- learning rate of the weight w
%          mu1 --- the Lagrange multiplier for constraint g1(y)
%          mu2 --- the Lagrange multiplier for constraint g2(w)
%       lambda --- the Lagrange multiplier for constraint h
%  OverValue --- stopping criterion, say, 0.0001. When the changed value 
%                of weighted vector is less than it, then the algorithm stops
%    maxIter --- maximum iterations for estimating each independent component
%
% See also:
%    genPulseRef     genRectangleRef     cICAdemo
%
% Reference:
%    [1] Wei Lu, Jagath C. Rajapakse: ICA with reference. ICA 2001
%    [2] Zhi-Lin Zhang, Morphologically Constrained ICA for Extracting Weak Temporally 
%        Correlated Signals, Neurocomputing 71(7-9) (2008) 1669-1679
%
% Author: Zhilin Zhang
%         z4zhang@ucsd.edu
%
% version: 1.0     Date: Dec.14,2008
% 


fprintf('Starting cICA for extracting the desired source signal ..\n');
[ICnum, IClen]=size(X);

w = w0;%initial w
oldw = w;

i = 1;
loop = 1;

meanref=ones(size(ref))*mean(ref);
stdref=std(ref);
refcankao=X*((ref-meanref)/stdref)'/IClen;

meanreft=ones(size(reft))*mean(reft);
stdreft=std(reft);
reftcankao=((reft-meanreft)/stdreft)/sqrt(ICnum);
 
thr1=zeros(200,1);thr2=zeros(200,1);
while (i == 1)
    
    learningRate =0.98^i;
            gamma = 0.1*4^(i-1);
            gamma1 =0.1*4^(i-1);
            gamma2 =0.2*4^(i-1);
           oldw=w;
    y = w'*X;       
    % calculate the first order deviation of the Lagarange function
    std_y = std(y);                          % standard deviation
    v_Gaus = normrnd(0, std_y, 1, IClen);    % Gaussian signal with the same mean and variance
    rou = mean( tanh(y) - tanh(v_Gaus));     %select nonlinear function tanh(x)  
    
     L1 = sign(rou) * ( X * tanh(y)')/IClen + mu1 * refcankao ...
        + mu2*reftcankao -lambda * ( X * y')/IClen;%select the variance as the referene
    
    % related to the second order deviation of the Lagarange function
    Sw = sign(rou) * mean(1-tanh(y).^2) - lambda;

    % compute the autocorrelation matrix Rxx
    Rxx = X(:,[1:end]) * X(:,[1:end])' / IClen;

    % update of the weight vector
    w = w - learningRate * inv(Rxx) * L1 / Sw; % 
    if i==1
        thr2(i,1)=threshold;
    elseif i==2
        thr2(i,1)=((reftt-meanreftt)'/stdreftt)*(w/sqrt(ICnum))+0.001;
    else
       if abs(thr2(i-1,1)-thr2(i-2,1))>0.01
           thr2(i,1)=((reftt-meanreftt)'/stdreftt)*(w/sqrt(ICnum))+0.001;
       elseif abs(thr2(i-1,1)-thr2(i-2,1))<=0.01
           thr2(i,1)=((reftt-meanreftt)'/stdreftt)*(w/sqrt(ICnum))-0.01;
       end
    end
    if i==1
        thr1(i,1)=threshold;
    elseif i==2
        thr1(i,1)=y*(ref-meanref)'/(stdref*IClen)+0.001;
    else
       if abs(thr1(i-1,1)-thr1(i-2,1))>0.01
           thr1(i,1)=y*(ref-meanref)'/(stdref*IClen)+0.01;
       elseif abs(thr1(i-1,1)-thr1(i-2,1))<=0.01
           thr1(i,1)=y*(ref-meanref)'/(stdref*IClen)-0.01;
       end
    end
    %thr = threshold;
    g1 = thr1(i,1)- y*(ref-meanref)'/(stdref*IClen);    % corresponds to the inequality constraint，mu的迭代关系
    mu1 = max(0, mu1 + gamma1 * g1);
%     
    g2 = thr2(i,1)-((reft-meanreft)'/stdreft)*(w/sqrt(ICnum)); 
    mu2 = max(0, mu2 + gamma2 * g2);
    % update of the parameter lambda
    h = mean(y.^2) - 1;                    % corresponds to the equality constraint lambda的迭代关系
    lambda = lambda + gamma * h;

    % decide whether the algorithm has converged or not
    wchange = norm(w-oldw);
    fprintf('No.%d iteration: change in w is %g\n',loop, wchange);
    if wchange < OverValue
        fprintf('Converged after %d iteration\n',loop);
        flag = 0;
        break
    end

%     if loop >= maxIter
%         fprintf('After %d iteration, still cannot convergent.\n',loop);
%         flag = 0;
%     end

    oldw = w;
    loop = loop + 1;

end

% output
IC = (w'* X)';

fprintf('End of cICA algorithm !\n');
