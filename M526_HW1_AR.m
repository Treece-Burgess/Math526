clear all;

load('streambed_data.mat')

u_func=@(z)(abs(z)<40).*(-15*exp(-(1/900)*z.^2))+(abs(z)>=40).*(z.*0); %mean(y)*z.^0;%(abs(z)<40).*(-15+0.01.*z.^2)+(abs(z)>=40).*(z.*0); %(abs(z)<40).*(-17*cos((z.*pi)./80))+(abs(z)>=40).*(0); 
lambda=1;
l=5;
n=101;
s=100000;                         
t_p=linspace(-100,100,n);

C_tp=lambda^2*exp(-(t_p'-t_p).^2/(2*l^2));
C_d=lambda^2*exp(-(x'-x).^2/(2*l^2));
C_dtp=lambda^2*exp(-(x'-t_p).^2/(2*l^2));
C_tpd=lambda^2*exp(-(t_p'-x).^2/(2*l^2));

L_pr=chol(C_tp);
u_pr=u_func(t_p);

prior=ones(s,1)*u_func(t_p)+normrnd(0,1,s,n)*L_pr;

% figure()
% hold on
% plot(t_p,prior, 'Color', "#888c89",'LineWidth',0.5)
% plot(t_p,mean(prior), 'Color','m','LineWidth',0.5)
% errorbar(x,y,d,'o','Color','black','LineWidth',1.5)
% title('Input Data vs Prior')
% ylabel('Output (y)')
% xlabel('Input (x)')
% hold off

G=inv(C_d+diag(d.^2));
U_pt=u_func(t_p)+(y-u_func(x))*G'*C_dtp;
C_pt=C_tp-C_tpd*G*C_dtp;
L_pt=chol(C_pt);
post=ones(s,1)*U_pt+normrnd(0,1,s,n)*L_pt;

% figure()
% hold on
% plot(t_p,post, 'Color', 'y', 'LineWidth', 0.5)
% plot(t_p,mean(post), 'Color', 'r', 'LineWidth', 0.5)
% errorbar(x,y,d,'o', 'Color', 'black', 'LineWidth', 1.5)
% title('Input Data vs Posterior')
% ylabel('Output (y)')
% xlabel('Input (x)')
% hold off
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


