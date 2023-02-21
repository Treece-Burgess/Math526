clear all;
load('streambed_data.mat')

%mean value
mean_func= @(z) exp(-(z.^2./500)).*(-30/sqrt(pi));

%parameters
lambda=0.25;
l=0.25;
%generate points 
test_pnts = linspace(-100,100,125);
%quadratic
q = @(x,y) lambda^2*exp(-(x-y).^2/(2*l^2));
%error
err = @(d) d.^2;
%Covariance matrix being passed
for i=1:length(test_pnts)
    for j=1:length(test_pnts)
        v(i,j)=q(test_pnts(i),test_pnts(j));
    end
end

%[U,S,∼] was the ~ symbol to just show it is an irrelevant value?
[U,S]=svd(v); 
L=U*sqrt(S);

%Drawing Monte Carlo Samples, how would you draw the samples?
for i=1:length(test_pnts)
    f(i,1)=normrnd(0,sqrt(1));
    mean_prior(i,1)=mean_func(test_pnts(i));

   %r(:,i) = mean + L * f;
end

r_prior=mean_prior + L * f;
%r_post=mean_post + L * f;


figure()
hold on
plot(test_pnts,r_prior)
errorbar(x,y,d,'o')
hold off

%Posterior
for i=1:length(test_pnts)
    for j=1:length(x)
        c_star_sharp(i,j)=q(test_pnts(i),x(j));
        c_sharp_sharp(i,j)=q(x(i),x(j));
        v_sharp_sharp(i,j)=err(d)*err(d)';
        c_sharp_star(i,j)=q(x(i),test_pnts(j));
    end
end

%Calculating Covariance Matrix Posterior
g_star_sharp=c_star_sharp*(c_sharp_sharp+v_sharp_sharp)^-1;
post_v= v - g_star_sharp*c_sharp_star;
[U_post,S_post,trash]=svd(post_v); 
L_post=U_post*sqrt(S_post);
%Calculating Posterior Mean
for i=1:length(test_pnts)
    f(i,1)=normrnd(0,1);
    mean_post(i,1)=q(x(i),x(j));
end

U=mean_prior+g_star_sharp*(y'-mean_post);

r_post=U + L_post*f;

e_prior=abs(y'-r_prior);
e_post=abs(y'-r_post);



figure()
hold on
%plot(test_pnts,r_prior,'-o','-r','LineWidth',2.0)
plot(test_pnts,r_prior,'-o','Color','red')
errorbar(x,y,d,'o')
plot(test_pnts,r_post,'-o','Color','green')
%plot(test_pnts,r_post,'-o','-g','LineWidth',2.0)
hold off