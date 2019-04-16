%this script investigates the weak convergence of the SGLD method 
%applied to a simple Gaussian model
clear all;clc

%creating the data
N=10^3; %this is how many data points we have
sdX=1; %the variance of our points
sdTheta=1; %the variance of the prior
sd=1; randn(sd);% Choose random number seed
theta=sdTheta*randn; %theta is N(0,sdTheta^2)
X=theta+sdX*randn(N,1); %these are the data points




%we will now define the mean and the variance of our posterior
muP=sum(X)/((sdX/sdTheta)^2+N);
sdP=1/sqrt(1/sdTheta^2+N/sdX^2);


%the SDE that we are solving is of the form
% dX=-0.5*(X-mu_p)/sigma_p^2 dt +dW

a=0.5/sdP^2; %this here is giving us the stiffness (increases a lot with N)

par=[0.1,0.25,0.5,0.75];

h=par/a;
n=10; %this how many points we will be taking into account in the calculation of the likelihood
vvar=N^2*(N-n)/(n*(N-1))*var(X);  %variance of the estimator
tic;
for l=1:length(h)
    T=0.02; %this is enough time to reach to equilibrium
    t=[0:h(l):T];
    M=5*10^4;
    K=length(t);
    
    
    
    Z=muP*ones(M,1);    %this is where  we are storing our information 
    sub=zeros(M,1);     %this is where  we are storing the resampling (probably a better way of doing this)
    Z1=muP*ones(M,1); 
    Z2=muP*ones(M,1);
    for i=1:length(t)
    
    
        for j=1:M %doing the resampling
            k=randsample(X,n); %without replacement
            sub(j)=N*mean(k);
        end
    
   
    Z=Z-h(l)*(0.5*Z/sdP^2-0.5*sub/sdX^2) +sqrt(h(l))*randn(M,1);
    %using all the points
    Z1=Z1-h(l)*(0.5*Z1/sdP^2-0.5*muP) +sqrt(h(l))*randn(M,1);
    Z2=Z2-h(l)*(0.5*Z2/sdP^2-0.5*sub/sdX^2) +sqrt(h(l)-h(l)^2/4*vvar)*randn(M,1);
    end
    %we know what the second moment should be for an OU at every time t
error(l)=abs(sdP^2*(1-exp(-2*T/sdP^2))-var(Z));
error_full(l)=abs(sdP^2*(1-exp(-2*T/sdP^2))-var(Z1));
error_full1(l)=abs(sdP^2*(1-exp(-2*T/sdP^2))-var(Z2));%msgld

end
toc;

%We are very close to equilibrium so in principle our error should be
%similar in size with the one we are calculating

%the variance of the estimator is (here we are without replacement)


error1=abs(sdP^2- (1+h*vvar)./(2*a-a^2*h));% sgld and the truth posterior
error2=abs(sdP^2- 1./(2*a-a^2*h)); % euler and the truth posterior
error3=abs(sdP^2- (1./(2*a-a^2*h)+h.*h*vvar.*vvar./(4*(2*a-a^2*h))))

% %figure showing error in SGLD being larger
% figure(1), loglog(h,error,'r--*',h,error_full,'b--o')
% title('figure showing error in SGLD being larger')
% %reality check when the full system is used between true and numerical
% %error
% figure(2), loglog(h,error_full,'b--o',h,error2,'r--*')
% title('euler method variance numerical and variance of estimator truth')
% 
% %reality check when the full system is used between true and numerical
% %error
% figure(3), loglog(h,error,'b--o',h,error1,'r--*')
% title('sgld method variance numerical and variance of estimator truth')
% 
%reality check when the full system is used between true and numerical
%error
% figure(4), loglog(h,error,'b--o',h,error_full1,'r--*')
% title('sgld,msgld')
% figure(5), loglog(h,error1,'r--*',h,error3,'b--o',h,error2,'k--^')
% title('n=10,the bias of sgld and msgld')
% xlabel('stepsize');ylabel('bias');
% legend('bias of sgld','bias of msgld','bias of Euler-Maruyama')
% 
% figure(6),loglog(h,error,h,error1,h,error2,h,error3)
% xlabel('stepsize');ylable('bias');
% legend('error','error1','error2','error3')

%%
%% the above is about the true one, the below is about our experiment 

format_str = {'Interpreter', 'latex','FontSize',30};
set(0, 'DefaultAxesFontSize',30);
set(0,'DefaultLineLineWidth', 4);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
figure(1),loglog(h,error,'-x',h,error_full,'-o',h,error_full1,'-^','MarkerSize',12);
xlabel('h', format_str{:});ylabel('Bias', format_str{:})
title('Experiment Bias of SGLD, mSGLD and Euler (n=10)', format_str{:})
legend('SGLD','Euler','mSGLD', format_str{:},'Location','NorthWest')
