x0 = [1.0,2.0,2.0,0.3,0.4,0.5];
global n sample
mu = [1,2,3]';
sig_matrix = diag([0.4^2,0.3^2,0.2^2]);
sample = mvnrnd(mu,sig_matrix,1000);
n = length(sample);
fun = @objfun;
[para_est,fval] = fminunc(fun,x0);

function s = objfun(para)
    global n sample
    mu_est = [para(1) para(2) para(3)];
    sig_est = diag([para(4).^2 para(5).^2 para(6).^2]);
    s = 0;
    for i = 1:n
       temp = (sample(i,:)-mu_est)*inv(sig_est)*(sample(i,:)-mu_est)';
       s = s+temp;
    end
    s = 1/2*s;
    s = s+n/2*log(det(sig_est));
    
end