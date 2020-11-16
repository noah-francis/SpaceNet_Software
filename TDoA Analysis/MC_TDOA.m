function [XYZ_mean,MAGstd_XYZ ] = MC_TDOA(k,error,Satellite,Rs,titleSTR)
%UNTITLED Summary of this function goes here
%   Monte Carlo simulation for TDoA uncertainty analysis with errors in
%k-number of simulations
%error- error in [err_timing, err_XYZpositions]
% True Position of the satellite in ECEF coordinates.

% Constants
c=299792458; %m/s

M=length(Rs);
XYZ=zeros(k,3);
err=zeros(k,3);
DOP=zeros(k,1);
M_mat=zeros(4,16,k);

titleSTR= strcat(titleSTR,'with \sigma_t =',num2str(error(1)*10^9),'ns and \sigma_r =', num2str(error(2)),'m');

for i=1:k
    err_t(i,:)=error(1)*randn(1,M);
    err_r(:,:,i)=error(2)*randn(size(Rs));
    [XYZ(i,:),err(i,:),DOP(i)]=TDOA(Satellite',Rs,c,err_r(:,:,i),err_t(i,:)',error(1));
end
Trash=find(abs(XYZ(:,1))>10^15);
err_t(Trash,:);
% err_r(Trash)
XYZ=XYZ(~isoutlier(XYZ(:,1)),:);
% % err(
XYZ=XYZ*10^(-3);

Satellite=Satellite*10^(-3);

% MAgnitude pop pop
MAG_XYZ=vecnorm(XYZ,2,2);
MAG_Sat=norm(Satellite);
% Mean
XYZ_mean=mean(XYZ);

MAGstd_XYZ=std(MAG_XYZ);

%Covariance Matrix

P=cov(XYZ);


sigma = sqrt( diag(P) );
% 
% Ecc=max(sigma);
% theta=linspace(0,2*pi,100);
% X_ell=2*sigma(1)*cos(theta);
% Y_ell=2*sigma(2)*sin(theta);
% Z_ell=2*sigma(3)*sin(theta);
% ANG_XY=atan2(XYZ_mean(2),XYZ_mean(1));
% 
% XY_plot=XYZ_mean(1)+X_ell*cos(ANG_XY)-Y_ell*sin(ANG_XY);
% YX_plot=XYZ_mean(2)+X_ell*sin(ANG_XY)+Y_ell*cos(ANG_XY);
% figure()
% plot(XY_plot,YX_plot)%2*XY_plot,2*YX_plot,3*XY_plot,3*YX_plot)
% hold on;
% plot(XYZ(:,1),XYZ(:,2),'.b')
figure;
subplot(2,2,1)
error_ellipse([XYZ(:,1),XYZ(:,2)])
grid on
ylabel('Y [km]')
subplot(2,2,3)
error_ellipse([XYZ(:,1),XYZ(:,3)])
grid on
xlabel('X[km]')
ylabel('Z [km]')
subplot(2,2,4)

error_ellipse([XYZ(:,2),XYZ(:,3)])
xlabel('Y [km]')


% Ecc=max(COV(i,:));
% theta=linspace(0,2*pi,100);
% POS=2*COV(i,1)*cos(theta);
% VEL=2*COV(i,2)*sin(theta);
% ANG=atan2(Epe(i,3),Epe(i,2));
% Pos=Epe(i,2)+POS*cos(ANG)-VEL*sin(ANG);
% Vel=Epe(i,3)+POS*sin(ANG)+VEL*cos(ANG);

%% Plot Histograms
figure()
sgtitle(titleSTR)
subplot(3,1,1)
histogram(XYZ(:,1),100);
hold on
grid on
xline(Satellite(1),'r','linewidth',2);
xline(XYZ_mean(1),'c','linewidth',2);
hold off
xlabel('TDoA X position Estimate [km]')
ylabel('Counts')
subplot(3,1,2)
histogram(XYZ(:,2),100);
hold on
grid on
xline(Satellite(2),'r','linewidth',2);
xline(XYZ_mean(2),'c','linewidth',2);
hold off
xlabel('TDoA Y position Estimate [km]')
ylabel('Counts')

subplot(3,1,3)
histogram(XYZ(:,3),100);
hold on
grid on
xline(Satellite(3),'r','linewidth',2);
xline(XYZ_mean(3),'c','linewidth',2);
hold off
xlabel('TDoA Z position Estimate [km]')
ylabel('Counts')
legend('Counts','True','Mean of TDoA outputs','location','best')


figure()

histogram(MAG_XYZ,100)
hold on
grid on
xline(MAG_Sat,'r','linewidth',2);
xline(mean(MAG_XYZ),'c','linewidth',2);
xline(mean(MAG_XYZ)-3*MAGstd_XYZ,'m','linewidth',2);
xline(mean(MAG_XYZ)+3*MAGstd_XYZ,'m','linewidth',2);
legend('Estimates','True Magnitude','Mean TDoA Estimate','3\sigma','location','Northeast')
xlabel('TDoA position Magnitude Estimate [km]')
ylabel('Counts')

title(titleSTR)

figure()

histogram((MAG_XYZ-MAG_Sat),100)
hold on
grid on
% xline(MAG_Sat,'r','linewidth',2)
% xline(mean(MAG_XYZ),'c','linewidth',2)
xline(-100,'m','linewidth',2);
xline(100,'m','linewidth',2);
legend('|R_e_s_t-R_T_r_u_e|','100 km Reguirment','location','Northwest')
xlabel('R_e_r_r_o_r [km]')
ylabel('Counts')
if 5*MAGstd_XYZ<100
    xlim([-150 150])
else
    xlim([-5*MAGstd_XYZ 5*MAGstd_XYZ])
end
title(titleSTR)


percentFail=(length(find(abs(MAG_XYZ-MAG_Sat)>=100))+k-length(MAG_XYZ))/k*100;

fprintf('%0.2f%% of the Position estimates were outside the 100 km requirement \n',percentFail)


end

function[XYZ,err,DOP]=TDOA(p_T,P,c,err_r,err_t,sig_t)
% This is based on TDoA file from MATLAB file exchange

method=1;


trials=100;
M=length(P);
in_est_error=2000;
% Induce Noise into Sensor Position SEED For repeatablility
% rng(40)


k_vec=[P(1:3,1)',P(1:3,2)',P(1:3,3)',P(1:3,4)']';

guess=P(1:3,1)+ P(1:3,1)/norm(P(1:3,1))*575e3;
for k=1:trials
%     err_r=sig_r*randn(size(P));
%     err_t=sig_t*randn(M-1,1);

%sensor positon vectors
% P = zeros(3,M); %includes all sensor positon vectors
% for ii=1:M
%     P(:,ii)=range_s*2*(rand(3,1)-0.5);
% end
%  
% p_T = range_T*2*(rand(3,1)-0.5);   %target positon vector
% %finding TOAs 
dummy = repmat(p_T,1,M)-P;
toa = zeros(M,1);   %includes all toa information
for ii = 1:M
    % Calculates time of arrivial and round to the nearest 420*10^(-7)
%     toa(ii) =norm(dummy(:,ii))/c-mod(norm(dummy(:,ii))/c,sig_t);    
    toa(ii) =norm(dummy(:,ii))/c/sig_t;
    toa(ii)=floor(toa(ii))*sig_t;
end
tdoa = toa-toa(1); %tdoa(1)=[];

P=P+err_r;
tdoa = tdoa +err_t;


% err_tMeters=[0;c*err_t];
Pt=(c*420*10^(-9))^2*eye(4);

    
%%% Taylor Series Expansion Solution
p_T_0 = guess + in_est_error*randn(3,1);    %initial estimate with some error (penalty term)
d = c*tdoa;
f = zeros(M,1);
del_f = zeros(M,3);

if method==1
    


for ii=1:M
   f(ii)=norm(p_T_0-P(:,ii))-norm(p_T_0-P(:,1)); 
   del_f(ii,1) = (p_T_0(1)-P(1,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(1)-P(1,1))*norm(p_T_0-P(:,1))^-1;
   del_f(ii,2) = (p_T_0(2)-P(2,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(2)-P(2,1))*norm(p_T_0-P(:,1))^-1;
   del_f(ii,3) = (p_T_0(3)-P(3,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(3)-P(3,1))*norm(p_T_0-P(:,1))^-1;    

end
%    del_f(1,:)=(p_T_0 - P(:,1))'*norm(p_T_0 - P(:,1))^(-1) - p_T_0'/norm(p_T_0);
   x_nonlin = pinv(del_f)*(d-f)+p_T_0;
   guess=x_nonlin;
   X=x_nonlin;
else
    A= -2 * [P(:,2)'- P(:,1)', d(2);
             P(:,3)'- P(:,1)', d(3); 
             P(:,4)'- P(:,1)', d(4)];
    b=[ d(2)^2-norm(P(:,2))^2+norm(P(:,1))^2;
        d(3)^2-norm(P(:,3))^2+norm(P(:,1))^2;
        d(4)^2-norm(P(:,4))^2+norm(P(:,1))^2];
%     Ainv=A'*A;
%     B=A'*b;
    X=pinv(A)*b;
%     X(1:3)=X(4)*X(1:3)/norm(X(1:3))+ X(1:3);
    
end
P=P-err_r;


% 
% rmse(k) = norm(p_T-X)^2;
end

   
% f=[0;f];
% df=gradient(f);
% dy=c*gradient([0;tdoa]);
% dk=gradient(k_vec);
% dr=gradient([0;p_T]);


% tdoa=[0;tdoa];

% dfdr=[0,0,0,0;
%        zeros(3,1),df./dr'];
% dfdy=[0,0,0,0;
%        zeros(3,1),df./dy'];
% dfdk=[zeros(1,12);
%        zeros(3,1),df./dk'];

%    dfdr=df./dr';
%    dfdy=df./dy';
%    dfdk=df./dk';

H=1/c*[ (p_T_0 - P(:,1))'*norm(p_T_0 - P(:,1))^(-1) - p_T_0'/norm(p_T_0);
    (p_T_0 - P(:,2))'*norm(p_T_0 - P(:,2))^(-1) - p_T_0'/norm(p_T_0);
    (p_T_0 - P(:,3))'*norm(p_T_0 - P(:,3))^(-1) - p_T_0'/norm(p_T_0);
    (p_T_0 - P(:,4))'*norm(p_T_0 - P(:,4))^(-1) - p_T_0'/norm(p_T_0)];

% M_mat=Covariace4TDoA(tdoa,k_vec,X,420*10^(-9));

DOP=sqrt(trace(pinv(H'*pinv(Pt)*H)));
XYZ=X;
err=p_T-X(1:3);
% if norm(err)>1e6
%   XYZ=NaN;
%   err=NaN;
%   DOP=NaN;
% end
end
