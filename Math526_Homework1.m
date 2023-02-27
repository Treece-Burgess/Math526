clear all;

%Section for arrays to store values
arr_neg20=zeros();
arr_neg10=zeros();
arr_zero=zeros();
arr_pos10=zeros();
arr_pos20=zeros();
pos=[51, 57, 63, 69, 75;];

%Loading in data for problem
load('streambed_data.mat')

%Defining anonymous function for the mean
mean_func=@(z) exp(-(z.^2./500)).*(-30/sqrt(pi));

%Definining lambda and l
lambda=0.25;
l=0.001;

%Generating Test Points
test_pnts=linspace(-100,100,125);

%Definining anonymous function for quadratic exponential
q=@(x,y) lambda^2*exp(-(x-y).^2/(2*l^2));

%Defining anonymous function for error 
err=@(d) d.^2;

%%Section below is for computing the prior

%Calculating covariance matrix for the prior
for i=1:length(test_pnts)
    for j=1:length(test_pnts)
        v(i,j)=q(test_pnts(i),test_pnts(j));
    end
end

%Performing SVD to obtain L for the prior
[U_prior,S]=svd(v); 
L=U_prior*sqrt(S);

%Calculating mean for the prior
for i=1:length(test_pnts)
    %f(i,1)=normrnd(0,sqrt(1));
    mean_prior(i,1)=mean_func(test_pnts(i));
end

figure()
n=0;
while n < 20

for i=1:length(test_pnts)
    f(i,1)=normrnd(0,sqrt(1));
end
r_prior=mean_prior + L*f;
hold on
txt = ['Monte Carlo Sample ',num2str(n)];
plot(test_pnts,r_prior,'DisplayName',txt,'LineWidth',2.0)
n = n+1;
end
errorbar(x,y,d,'o','Color','black','DisplayName','Input Data','LineWidth',2.0)
title('Input Data vs Prior')
ylabel('Output (y)')
xlabel('Input (x)')
legend('Location','east')
hold off

%%Section below is for computing the posterior

%Computing matrices for the posterior calculation
for i=1:length(test_pnts)
    for j=1:length(x)
        c_star_sharp(i,j)=q(test_pnts(i),x(j));
        c_sharp_sharp(i,j)=q(x(i),x(j));
        v_sharp_sharp(i,j)=err(d)*err(d)';
        c_sharp_star(i,j)=q(x(i),test_pnts(j));
    end
end

%Computing the covariance matrix 
g_star_sharp=c_star_sharp*(c_sharp_sharp+v_sharp_sharp)^-1;
post_v= v-g_star_sharp*c_sharp_star;

%Performing SVD to obtain L for the posterior
[U_post,S_post]=svd(post_v); 
L_post=U_post*sqrt(S_post);

%Calculating mean for the posterior
for i=1:length(test_pnts)
    %f(i,1)=normrnd(0,1);
    mean_post(i,1)=mean_func(x(i));
end

mu=mean_prior+g_star_sharp*(y'-mean_post);

figure()
n=0;
while n < 20
for i=1:length(test_pnts)
    f(i,1)=normrnd(0,sqrt(1));
end
r_post=mu + L_post*f;
hold on
txt = ['Monte Carlo Sample ',num2str(n)];
plot(test_pnts,r_post,'DisplayName',txt,'LineWidth',2.0)
n = n+1;
end

errorbar(x,y,d,'o','Color','black','DisplayName','Input Data','LineWidth',2.0)
title('Input Data vs Posterior')
ylabel('Output (y)')
xlabel('Input (x)')
legend('Location','east')
hold off

%%Number 4
%Analytic
x_let=[-20,-10,0,10,20;];
num_4=[51, 57, 63, 69, 75;];
for i=1:5
cdf4(i)=0.5*(1-erf((-10-mu(num_4(i)))/(sqrt(2)*sqrt((L_post(num_4(i),num_4(i)))^2))));
end
%Monte Carlo Sampling
n=1;
while n < 10001
for i=1:length(test_pnts)
    f(i,1)=normrnd(0,sqrt(1));
end
r_post=mu + L_post*f;
arr_neg20(1,n) = r_post(51);
arr_neg10(1,n) = r_post(57);
arr_zero(1,n) = r_post(63);
arr_pos10(1,n) = r_post(69);
arr_pos20(1,n) = r_post(75);
n = n + 1;
end

1/10000 * sum(arr_neg20 > -10)
1/10000 * sum(arr_neg10 > -10)
1/10000 * sum(arr_zero > -10)
1/10000 * sum(arr_pos10 > -10)
1/10000 * sum(arr_pos20 > -10)

%%Number 5
arr_simul=zeros();
for i=1:10000
    if(arr_neg20(i) > -10 && arr_neg10(i) > -10 && arr_zero(i) > -10 && arr_pos10(i) > -10 && arr_post20(i) > -10)
        arr_simul(i) = 1;
    else
        arr_simul(i) = 0;
    end
end

1/10000 * sum(arr_simul)
