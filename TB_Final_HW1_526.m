clear all;

load('streambed_data.mat')

u_func=@(z) (abs(z)<40).*(-15*exp(-(1/1000)*z.^2))+(abs(z)>=40).*(z.*0);%(abs(z)<40).*(-15+0.01.*z.^2)+(abs(z)>=40).*(z.*0);%(abs(z)<40).*(-15*exp(-(1/1000)*z.^2))+(abs(z)>=40).*(z.*0); %mean(y)*z.^0;% %(abs(z)<40).*(-17*cos((z.*pi)./80))+(abs(z)>=40).*(0); 
lambda=5;
l=5;
n=101;
s=10000;                         
t_p=linspace(-100,100,n);

C_tptp=lambda^2*exp(-(t_p'-t_p).^2/(2*l^2));
C_dd=lambda^2*exp(-(x'-x).^2/(2*l^2));
C_dtp=lambda^2*exp(-(x'-t_p).^2/(2*l^2));
C_tpd=lambda^2*exp(-(t_p'-x).^2/(2*l^2));

[V,S,~]=svd(C_tptp);
%L_pr=(V*S^0.5)';
L_pr=chol(C_tptp); 
%Cholesky for 4 & 5 since the matrix is PD
u_pr=u_func(t_p);
prior=ones(s,1)*u_func(t_p)+normrnd(0,1,s,n)*L_pr;
% 
figure()
hold on
h1=plot(t_p,prior, 'Color', "#00FFFF",'LineWidth',0.5);
h2=plot(t_p,mean(prior), 'Color','#FF0000','LineWidth',2.0);
h3=errorbar(x,y,d,'o','Color','black','LineWidth',1.5);
title('Input Data vs Prior (lambda=1,l=5)')
ylabel('Output (y)')
xlabel('Input (x)')
legend([h1(1),h2(1),h3(1)],'Prior','Prior MAP Estimate','Input Data','Location','southeast')
hold off

G=inv(C_dd+diag(d.^2));
U_pt=u_func(t_p)+(y-u_func(x))*G'*C_dtp;
C_pt=C_tptp-C_tpd*G*C_dtp;
[V,S,~]=svd(C_pt);
%L_pt=(V*sqrt(S))';
L_pt=chol(C_pt); %Cholesky for 4&5 when the matrix is PD
post=ones(s,1)*U_pt+normrnd(0,1,s,n)*L_pt;

figure()
hold on
h1=plot(t_p,post, 'Color', '#0000FF', 'LineWidth', 1.0);
h2=plot(t_p,mean(post), 'Color', 'r', 'LineWidth', 2.0);
h3=errorbar(x,y,d,'o', 'Color', 'black', 'LineWidth', 1.5);
title('Input Data vs Posterior (lambda=5,l=1)')
ylabel('Output (y)')
xlabel('Input (x)')
legend([h1(1),h2(1),h3(1)],'Posterior','Posterior MAP Estimate','Input Data','Location','southeast')
hold off
% 
figure()
hold on
plot(t_p,prior, 'Color', "#888c89",'LineWidth',0.25)
plot(t_p,mean(prior), 'Color','m','LineWidth',0.25)
plot(t_p,post, 'Color', 'y', 'LineWidth', 0.25)
plot(t_p,mean(post), 'Color', 'r', 'LineWidth', 0.25)
errorbar(x,y,d,'o','Color','black','LineWidth',1.5)
title('Input Data with Prior and Posterior')
ylabel('Output (y)')
xlabel('Input (x)')
hold off

% Part 4 and 5
a = find(t_p == -20);
b = find(t_p == -10);
c = find(t_p == 0);
d = find(t_p == 10);
e = find(t_p == 20);
alpha = [a, b, c, d, e;];

x_alpha=[-20, -10, 0, 10, 20];
u_alpha=[U_pt(a),U_pt(b),U_pt(c),U_pt(d),U_pt(e)];
sigma_D=[L_pt(a,a),L_pt(b,b),L_pt(c,c),L_pt(d,d),L_pt(e,e)];
sigma_alpha=diag(sigma_D);

for i=1:length(alpha)
    cdf4(i)=1-normcdf(-10,u_alpha(i),sigma_alpha(i,i));
    cdf_4(i)=0.5*(1-erf((-10-u_alpha(i))/(sqrt(2)*(sigma_alpha(i,i)))));
    
end

I_4=zeros(length(alpha),1);
for i=1:length(alpha)
    mc_4=zeros(s,1);
    for j=1:s
        if post(j,alpha(i)) > -10
            mc_4(j)=1;
        else            
        end
    end   
I_4(i)=(1/s)*sum(mc_4); 
end

cdf_5=mvncdf([-10,-10,-10,-10,- 10], u_alpha, sqrt(sigma_alpha.^2));
for i=1:s
    mc_5=zeros(s,1);
    
    if (post(i,a) > -10) && (post(i,b) > -10) && (post(i,c) > -10)...
            && (post(i,d) > -10) && (post(i,e) > -10)
        mc_5(i)=1;
    else
    end  
end
I_5=(1/s)*sum(mc_5);