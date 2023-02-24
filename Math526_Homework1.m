clear all;

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
[mu,S]=svd(v); 
L=mu*sqrt(S);

%Calculating mean for the prior
for i=1:length(test_pnts)
    f(i,1)=normrnd(0,sqrt(1));
    mean_prior(i,1)=mean_func(test_pnts(i));
end

%Monte Carlo Sample for the prior
r_prior=mean_prior + L * f;

%Generates figure for the input data and prior
figure()
hold on
title('Input Data vs Prior')
ylabel('Output (y)')
xlabel('Input (x)')
errorbar(x,y,d,'o','Color','red')
plot(test_pnts,r_prior,'-o','Color','blue')
legend('Input Data','Prior')
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
    f(i,1)=normrnd(0,1);
    mean_post(i,1)=q(x(i),x(j));
end

mu=mean_prior+g_star_sharp*(y'-mean_post);

%Monte Carlo Samples for the posterior
r_post=mu + L_post*f;


%Generates figure for the input data, prior, and posterior
figure()
hold on
title('Input Data vs Prior vs Posterior')
ylabel('Output (y)')
xlabel('Input (x)')
errorbar(x,y,d,'o','Color','red')
plot(test_pnts,r_prior,'-o','Color','blue')
plot(test_pnts,r_post,'-o','Color','green')
legend('Input data','Prior','Posterior')
hold off