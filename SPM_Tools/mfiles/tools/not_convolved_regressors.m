% Getting the convoluted vs non convoluted regresors, nonconvoluted
% regressors are stored in 

convolved=SPM.xX.X(:,2);
notconvolved=SPM.Sess.U.u;
figure;plot(notconvolved(:,2));title('NOT CONVOLVED SPM.Sess.U.u(:,2)')
figure;plot(convolved);title('CONVOLVED SPM.xX.X(:,2)')