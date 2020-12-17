%% Housekeeping
clc;clear;close all;
%% TDoA Code rev1.0
%{
By Keith Poletti 10/4/2020
This code aims to take the input time of arrivals, location of 4 sensor
units, and uncertainty in time sync/position. It will output a prediction
for the position of a satellite. 

This code is largely based on "UNCUED SATELLITE INITIAL ORBIT DETERMINATION USING
SIGNALS OF OPPORTUNITY" By Johnny L. Worthy III, Marcus J. Holzinger.

It will also use there given conditions as test cases.
%}
%%
c=299792458; %m/s
% define error in the time 
sig_t=500*10^(-9); % error in sec
sig_r=10; % error in position m

LineOfSite=1500*10^3;

% Access= xlsread('CSIM_To_All_SpaceNet.csv');
% Range=xlsread('CSIM_Station_Range.csv');
load CSIM_Position&Veloctiy_ECIJ2000.mat


S_t=cell2mat(finalCell);
% S_t=xlsread('CSIM_ECEF_Position.csv');

dates=S_t(:,1:6);
%gets position and converts to m from km
S_t=S_t(:,7:end)*10^3;
Xtrue_eci=S_t(:,1);
Ytrue_eci=S_t(:,2);
Ztrue_eci=S_t(:,3);

% Uncomment if data is new and there is no ECEF file
% lla = eci2lla(S_t(:,1:3),dates);
% RSatECEF=lla2ecef(lla);
% S_t_Ecef=RSatECEF;
% save('ECEF.mat','RSatECEF')
load ECEF.mat
S_t_Ecef=RSatECEF;

% S_t_Ecef=S_t(:,2:4)*10^3;

% Get sensor locations on earth center fixed and adds error
Rs=[lla2ecef([40.00888, -105.24774,1612]);
    lla2ecef([40.935583, -105.380917,1868]);
	lla2ecef([39.95121, -106.34978,2331]);
% 	lla2ecef([39.35773, -104.58592,1989]);
    lla2ecef([38.24047, -104.57511,2000])
%     lla2ecef([39.629142, -104.796563,1989])
    ];


%% finds relative distance for each sensor
TMP=vecnorm(S_t_Ecef-Rs(1,:),2,2);
r1=S_t_Ecef(abs(TMP)<LineOfSite,1:3);

TMP=vecnorm(S_t_Ecef-Rs(2,:),2,2);
r2=S_t_Ecef(abs(TMP)<LineOfSite,1:3);

TMP=vecnorm(S_t_Ecef-Rs(3,:),2,2);
r3=S_t_Ecef(abs(TMP)<LineOfSite,1:3);


TMP=vecnorm(S_t_Ecef-Rs(4,:),2,2);
r4=S_t_Ecef(abs(TMP)<LineOfSite,1:3);

[~,short_boi]=min([length(r1(:,1)),length(r2(:,1)),length(r3(:,1)),length(r4(:,1))]);



TMP=vecnorm(S_t_Ecef-Rs(short_boi,:),2,2);

% 
dates=dates(abs(TMP)<LineOfSite,:);
R=S_t(abs(TMP)<LineOfSite,:);


Xtrue=R(:,1);
Ytrue=R(:,2);
Ztrue=R(:,3);

R_ECEF=S_t_Ecef(abs(TMP)<LineOfSite,:);
ind_ECI=linspace(1,length(S_t),length(S_t))';
ind_ECI=ind_ECI(abs(TMP)<LineOfSite);
TMP=TMP(abs(TMP)<LineOfSite);

%% MONTE CARLO GONNA WIN All the money
k=10000;
fprintf('Worst Case: sig_t= %f ns sig_r= %f m\n',sig_t*10^9,sig_r)
%Finds the indices of different pass starts and stops
diffIND=find(abs(diff(ind_ECI))>1)+1;
% Find the Max range possible  given the line of sight and it's index
[~,WorstCase]=max(TMP);


% This worst cases will 100% happen at a beginning or end of a fly by so we
% must account for both cases
conditional=find(WorstCase==diffIND);
if ~isnan(conditional)
    %In this case the Worst case is at a start of a fly by
    start=WorstCase;
    stop=diffIND(conditional+1)-1;
   [~,mid]=min(TMP(start:stop));
   mid=mid+start;
    
else
    %In this case the Worst case is at a end of a fly by
    [~,start]=min(abs(WorstCase-diffIND));
    start=diffIND(start-1);
    stop=WorstCase;
    [~,mid]=min(TMP(start:stop));
    mid=mid+start;
end
[XYZHorizon1_Mean,MAGHorizon1_std]=MC_TDOA(k,[sig_t,sig_r],R_ECEF(start,:),Rs','First Horizon TDoA Monte Carlo');
[XYZClose_Mean,MAGClose_std]=MC_TDOA(k,[sig_t,sig_r],R_ECEF(mid,:),Rs','Closest to LASP TDoA Monte Carlo');
[XYZHorizon2_Mean,MAGHorizon2_std]=MC_TDOA(k,[sig_t,sig_r],R_ECEF(stop,:),Rs','Second Horizon TDoA Monte Carlo');

XYZ=[XYZHorizon1_Mean;
    XYZClose_Mean;
    XYZHorizon2_Mean];
%% plot Worst Case

figure
earth_sphere(50,'m')
hold on;
plot3(Rs(:,1), Rs(:,2),Rs(:,3),'m.','MarkerSize',7.75); hold on;
% plot3(shortstak(:,1),shortstak(:,2),shortstak(:,3),'r*');
% plot3(R_ECEF(1:120,1), R_ECEF(1:120,2),R_ECEF(1:120,3),'r*','MarkerIndices',1:5:120);
plot3(XYZ(:,1)*10^(3),XYZ(:,2)*10^(3),XYZ(:,3)*10^(3),'ms','MarkerSize',7.75);
% plot3(X, Y,Z,'ms','MarkerSize',7.75);
plot3(R_ECEF(start:stop,1),R_ECEF(start:stop,2),R_ECEF(start:stop,3),'c','linewidth',2)
% quiver3(Rs(1,1),Rs(1,2),Rs(1,3),TMP(1330,1),TMP(1330,2),TMP(1330,3))

legend('Earth','Sensor Positions', 'Target Estimation','Orbit')
grid on; 



%% Best Case 
% close all
fprintf('Best Case: sig_t= %f ns sig_r= %f m\n',sig_t*10^9,sig_r)
[~,BestCase]=min(TMP);

conditional=find(BestCase<=diffIND,1);
if isempty(conditional )
    start=diffIND(end);
    stop=length(TMP);
    mid=BestCase;
end

[XYZHorizon1_Mean_Best,MAGHorizon1_std]=MC_TDOA(k,[sig_t,sig_r],R_ECEF(start,:),Rs','First Horizon TDoA Monte Carlo');
[XYZClose_Mean_Best,MAGClose_std]=MC_TDOA(k,[sig_t,sig_r],R_ECEF(mid,:),Rs','Closest to LASP TDoA Monte Carlo');
[XYZHorizon2_Mean_Best,MAGHorizon2_std]=MC_TDOA(k,[sig_t,sig_r],R_ECEF(stop,:),Rs','Second Horizon TDoA Monte Carlo');

XYZ_Best=[XYZHorizon1_Mean_Best;
    XYZClose_Mean_Best;
    XYZHorizon2_Mean_Best];
% plot Best Case
figure
earth_sphere(50,'m')
hold on;
plot3(Rs(:,1), Rs(:,2),Rs(:,3),'m.','MarkerSize',7.75); hold on;
plot3(XYZ_Best(:,1)*10^(3),XYZ_Best(:,2)*10^(3),XYZ_Best(:,3)*10^(3),'ms','MarkerSize',7.75);

plot3(R_ECEF(start:stop,1),R_ECEF(start:stop,2),R_ECEF(start:stop,3),'c','linewidth',2)

legend('Earth','Sensor Positions', 'Target Estimation','Orbit')
grid on; 


close all
%% Residual vs flyby
for i=start:stop
    [~,err,~,~]=TDOA(R_ECEF(i,:)',Rs',c,sig_r,sig_t);
    
    D2LASP(i-start+1)=norm(R_ECEF(i,:)-Rs(1,:));
    residual(i-start+1)=norm(err)*10^-3;
    
end
figure();

plot(1:length(residual),residual)


% 
% 
% 
%% Test Based on Holzinger Paper
% X=[-0.9121 0.3707 0.6046]*6378*10^3;
% Senr=[lla2eci([33.7774, -84.3989,320],[2020,11,5,15,35,0]);
%     lla2eci([34.1374, -118.1256,250],[2020,11,5,15,35,0]);
% 	lla2eci([42.3608, -71.0928,30],[2020,11,5,15,35,0]);
% 	lla2eci([30.2852, -97.7348,158],[2020,12,5,15,35,0])];
% 
% 
%  [posHolz,errHolz,DOPHolz,MHolz]=TDOA(X',Senr',c,1,5e-9)
% COVmat=[eye(12),  zeros(12,4);
%          zeros(4,12),    (5e-9)^2*eye(4)];
%  
% [~,ind]=min(DOPHolz);
% PR=MHolz(:,:,ind)*COVmat*MHolz(:,:,ind)'
% Sig_XYZt_KM=sqrt(diag(PR))
% RMSERR=norm(sqrt(diag(PR)))
% DOPmin=DOPHolz
% 
%  
% figure
% earth_sphere(50,'m')
% hold on;
% plot3(Senr(:,1), Senr(:,2),Senr(:,3),'m.','MarkerSize',7.75); hold on;
% % plot3(shortstak(:,1),shortstak(:,2),shortstak(:,3),'r*');
% plot3(X(1), X(2),X(3),'r*','MarkerIndices',1:5:120);
% plot3(posHolz(1),posHolz(2),posHolz(3),'ms','MarkerSize',7.75,'MarkerIndices',1:5:120)
% % plot3(X, Y,Z,'ms','MarkerSize',7.75);
% %plot3(S_t_Ecef(1:5000,1),S_t_Ecef(1:5000,2),S_t_Ecef(1:5000,3),'c','linewidth',2)
% % quiver3(Rs(1,1),Rs(1,2),Rs(1,3),TMP(1330,1),TMP(1330,2),TMP(1330,3))
% 
% legend('Earth','Sensor Positions', 'Target Position', 'Target Estimation','Orbit')
% grid on; 

%%
% sig_t=linspace(30*10^(-9),420*10^(-9),150);
% 
% 
% for i=1:length(sig_t)
%    % Find master receiver
%     [pos1(i,:),err1(i,:),DOP1(i,:),M1(:,:,i)]=TDOA(R_ECEF(1,:)',Rs',c,sig_r,sig_t(i));
%     [pos2(i,:),err2(i,:),DOP2(i,:),M2(:,:,i)]=TDOA(R_ECEF(17,:)',Rs',c,sig_r,sig_t(i));
%     err_r1(i)=norm(err1(i,:))*10^(-3);
%     err_r2(i)=norm(err2(i,:))*10^(-3);
%     r(i)=norm(R_ECEF(i,:)-Rs(1,:))*10^(-3);
% end
% 
% MeanERR=mean(err_r2)
% err_STD=std(err_r2)
% figure
% plot(sig_t,err_r1)
% yline(100,'r','linewidth',2)
% grid on
% title('Horizon Position Error Estimate vs Timing Accuracy')
% xlabel('Timing Uncertainty [s]')
% ylim([0 250])
% ylabel('Error Vector Magnitude[km]')
% figure
% plot(sig_t,err_r2)
% yline(100,'r','linewidth',2)
% grid on
% title('Closest Point Position Error Estimate vs Timing Accuracy')
% xlabel('Timing Uncertainty [s]')
% ylim([0 250])
% ylabel('Error Vector Magnitude [km]')
% sig_t=420*10^(-9); % error in sec


clear err
%%
COVmat=[sig_r^2*eye(12),  zeros(12,4);
        zeros(4,12),      sig_t^2*eye(4)];
    
% Senr=[lla2ecef([33.7774, -84.3989,320]);
%       lla2ecef([34.1374, -118.1256,250]);
% 	  lla2ecef([42.3608, -71.0928,30]);
% 	  lla2ecef([30.2852, -97.7348,158])];
% sig_t=420*10^(-9); 

 for j=1:length(R(:,1))       
     
    %Find time from satellite to sensor unit
    % tau=All_access(:,3)/c;
%     Rs=[lla2eci([40.00888, -105.24774,1612],dates(j,:));
%     lla2eci([40.99764, -104.93957,1868],dates(j,:));
% 	lla2eci([39.95121, -106.34978,2331],dates(j,:));
% 	lla2eci([39.35773, -104.58592,1989],dates(j,:))];
    
    % Find master receiver
    [pos(j,:),err(j,:),DOP(j,:),M(:,:,j)]=TDOA(R_ECEF(j,:)',Rs',c,sig_r,sig_t);
    err_r(j)=norm(err(j,:))*10^(-3);
    r(j)=norm(R_ECEF(j,:)-Rs(1,:))*10^(-3);
    
    
    dcm = dcmeci2ecef('IAU-2000/2006',dates(j,:));
    Sen_ECI(:,:,j)=dcm\Rs';
    pos_ECI(j,:)=dcm\pos(j,:)';
    PR(:,:,j)=M(:,:,j)*COVmat*M(:,:,j)';
%     PR(1,1,j)=sig_t*c*DOP(j,:).^2;
    Sig_XYZt(j,:)=sqrt(diag(PR(:,:,j)));
    RMSERR(j,:)=norm(sqrt(diag(PR(:,:,j))))*10^(-3);


 end
%%
TRASH=find(err_r>10^5);
pos(TRASH,:)=[];
err_r(TRASH)=[];
DOP(TRASH)=[];
M(:,:,TRASH)=[];
r(TRASH)=[];
dates(TRASH,:)=[];
Sig_XYZt(TRASH,:)=[];
R(TRASH,:)=[];
ind_ECI(TRASH)=[];
pos_ECI(TRASH,:)=[];
Sen_ECI(:,:,TRASH)=[];

save('ECI_Sensor_Locations.mat','Sen_ECI')
% 
save('Covariance_of_Position.mat','PR')
% 
% 

save('TDoA_DOP.mat','DOP')

MeanERR=mean(err_r)
err_STD=std(err_r)
err_r=err_r';
ECI_Tdoa=[dates,pos_ECI,Sig_XYZt,err_r];
save('TDoA_ECI_out.mat','ECI_Tdoa')
% %     pos(:,1:3)=pos(:,1:3)+Rs(1,:);

COVmat=[sig_r^2*eye(12),  zeros(12,4);
         zeros(4,12),    sig_t^2*eye(4)];
ECI_True=[dates,R,ind_ECI];
save('ECI_True.mat','ECI_True')
%%
[~,ind]=min(DOP);
PR=M(:,:,ind)*COVmat*M(:,:,ind)'
Sig_XYZt_KM=sqrt(diag(PR))*10^(-3)
RMSERR=norm(sqrt(diag(PR)))*10^(-3)
DOPmin=DOP(ind)


% Position uncertainty on based TDoA in meter
OVERALL_TDOA_STD=DOP*c*sig_t*10^(-3);

figure;

plot(r,OVERALL_TDOA_STD,'.')
grid on
ylabel('c\sigma _t DOP [km]')
xlabel('Satellite Distance from LASP [km]')
%%
%{
To continue, Monte carlo simulations of different sensor locations in CO
    Be able to take multiple observations at different times in orbit
    See how increasing sensor units affects error.
%}

%%

% d1=vecnorm(R-Rs(1,:),2,2)*10^(-3);
% logic=d1>2768.66;
% d2=vecnorm(R-Rs(2,:),2,2)*10^(-3);
% logic=d2>2768.66+logic;
% d3=vecnorm(R-Rs(3,:),2,2)*10^(-3);
% logic=d3>2768.66+logic;
% d4=vecnorm(R-Rs(4,:),2,2)*10^(-3);
% logic=d4>2768.66+logic;
%%
% Rs=[lla2ecef([40.00888, -105.24774,1612]);
%     lla2ecef([40.99764, -104.93957,1868]);
% 	lla2ecef([39.95121, -106.34978,2331]);
% 	lla2ecef([39.35773, -104.58592,1989])];
figure
earth_sphere(50,'m')
hold on;
plot3(Rs(:,1), Rs(:,2),Rs(:,3),'m.','MarkerSize',7.75); hold on;
% plot3(shortstak(:,1),shortstak(:,2),shortstak(:,3),'r*');
plot3(R_ECEF(1:120,1), R_ECEF(1:120,2),R_ECEF(1:120,3),'r*','MarkerIndices',1:5:120);
plot3(pos(:,1),pos(:,2),pos(:,3),'ms','MarkerSize',7.75,'MarkerIndices',1:5:120);
% plot3(X, Y,Z,'ms','MarkerSize',7.75);
plot3(S_t_Ecef(1:5000,1),S_t_Ecef(1:5000,2),S_t_Ecef(1:5000,3),'c','linewidth',2)
% quiver3(Rs(1,1),Rs(1,2),Rs(1,3),TMP(1330,1),TMP(1330,2),TMP(1330,3))

legend('Earth','Sensor Positions', 'Target Position', 'Target Estimation','Orbit')
grid on; 
%%
figure
plot(r,err_r,'.')
yline(100,'r','linewidth',2);
grid on
title('Position Error Estimate vs Orbital Radius')
xlabel('Satellite Distance from LASP [km]')
ylim([0 500])
ylabel('Normalized Error Vector [km]')

%%
% clear r_sen err_sen
% R=lla2ecef([38.254762, -104.634134,575000]);
% sig_t=30*10^(-9);
% dAngle=linspace(0,12,100);
% k=1;
% Leg_str={'Earth', 'Target Position'};
% figure;
% earth_sphere(50,'m')
% hold on;
% plot3(R(1), R(2),R(3),'r*');
% for i=2:100
%     Rs=[lla2ecef([40.001262+dAngle(i), -105.269597,1661]);
%         lla2ecef([40.001262, -105.269597+dAngle(i),1661]);
%         lla2ecef([40.001262-dAngle(i), -105.269597,1661]);
%         lla2ecef([40.001262, -105.269597-dAngle(i),1661])];
%     [X,Y,Z,err_sen(i,:)]=TDOA(R,Rs,c,sig_r,sig_t);
%     pos_sen(i,:)=[X,Y,Z];
%     tmp=norm(Rs(1,:)-lla2ecef([40.001262, -105.269597,1661]));
%     r_sen(i,1)=tmp;
%     tmp=norm(Rs(2,:)-lla2ecef([40.001262, -105.269597,1661]));
%     r_sen(i,2)=tmp;
%     tmp=norm(Rs(3,:)-lla2ecef([40.001262, -105.269597,1661]));
%     r_sen(i,3)=tmp;
%     tmp=norm(Rs(4,:)-lla2ecef([40.001262, -105.269597,1661]));
%     r_sen(i,4)=tmp;
%     if mod(i,25)==0|| i==2
%         plot3(X, Y,Z,'ms','MarkerSize',7.75);
%         plot3(Rs(:,1), Rs(:,2),Rs(:,3),'.','MarkerSize',7.75); hold on;
%         Leg_str={Leg_str,strcat('r_x=', sprintf('%.6f',r_sen(i,1))), strcat('r_y=', sprintf('%.6f',r_sen(i,2)))};
%     end
% end
% % legend(Leg_str)
% r_sen=r_sen*10^(-3);
% %%
% figure;
% semilogx(r_sen(2:end,1),err_sen(2:end,2))
% yline(10,'r','linewidth',2)
% grid on
% title('Position Error Estimate vs Time accuracy')
% xlabel('Distance between LASP and Sensors [km]')
% ylabel('Normalized Error Vector [km]')
% ylim([0,50])
%%
function[SN1,SN2,SN3,SN4]=parser(Range,Access)
Access=Access(~isnan(Access(:,1)),:);
First=[find(Access(:,1)==1);length(Access(:,1))];
ind=cell(length(First)-1,1);
for i=1:length(First)-1
    v=First(i):First(i+1)-1;
    for j=1:length(v)
        tmp=find(Range(:,1)>=Access(v(j),2) & Range(:,1)<=Access(v(j),3));
        ind{i}=[ind{i};tmp];
    end
end
find(ind{1}==ind{2}==ind{3}==ind{4})
%define Each Radio
SN1=[ind{1},Range(ind{1},1),Range(ind{1},2)];
SN2=[ind{2},Range(ind{2},1),Range(ind{2},3)];
SN3=[ind{3},Range(ind{3},1),Range(ind{3},4)];
SN4=[ind{4},Range(ind{4},1),Range(ind{4},5)];


end

function[XYZ,err,DOP,M_mat]=TDOA(p_T,P,c,sig_r,sig_t)
% This is based on TDoA file from MATLAB file exchange

method=1;


trials=100;
M=length(P);
in_est_error=2000;
% Induce Noise into Sensor Position SEED For repeatablility
% rng(40)
err_r=sig_r*randn(size(P));
err_t=sig_t*randn(M,1);

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
    toa(ii) =(norm(dummy(:,ii))/c + err_t(ii))/sig_t;
    toa(ii)=floor(toa(ii))*sig_t;   
end
tdoa = toa-toa(1); %tdoa(1)=[];

Perr=P+err_r;
% tdoa = tdoa +err_t;


err_tMeters=[0;c*err_t];
Pt=(sig_t)^2*eye(4);

    
%%% Taylor Series Expansion Solution
p_T_0 = guess + in_est_error*randn(3,1);    %initial estimate with some error (penalty term)
d = c*tdoa;
f = zeros(M,1);
del_f = zeros(M,3);

if method==1
    


for ii=1:M
   f(ii)=norm(p_T_0-Perr(:,ii))-norm(p_T_0-Perr(:,1)); 
   del_f(ii,1) = (p_T_0(1)-Perr(1,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(1)-Perr(1,1))*norm(p_T_0-Perr(:,1))^-1;
   del_f(ii,2) = (p_T_0(2)-Perr(2,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(2)-Perr(2,1))*norm(p_T_0-Perr(:,1))^-1;
   del_f(ii,3) = (p_T_0(3)-Perr(3,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(3)-Perr(3,1))*norm(p_T_0-Perr(:,1))^-1;    

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

M_mat=Covariace4TDoA(tdoa,k_vec,X,sig_t);

DOP=sqrt(trace(pinv(H'*pinv(Pt)*H)));
XYZ=X;
err=p_T-X(1:3);
% if norm(err)>1e6
%   XYZ=NaN;
%   err=NaN;
%   DOP=NaN;
% end
end
% function[X,Y,Z,err]=TDOA(R,Rs,c,sig_r,sig_t)
% Rs_rand=Rs+sig_r*(-1 + 2.*randn(4,3));
% for i=1:length(Rs)
%     tau(i)=round(1/c*norm(Rs(i,:)-R),7);
% end
% err_t=[0;sig_t*randn(length(Rs)-1,1)];
% delT=tau-tau(1)+err_t';%[0;sig_t(j)*ones(length(Rs)-1,1)]';
% 
% 
% x=linspace( 2000*10^3+min(Rs_rand(:,1)) , max(Rs_rand(:,1))+2000*10^3,100);
% y=linspace( 2000*10^3+min(Rs_rand(:,2)) , max(Rs_rand(:,2))+2000*10^3,100);
% z=linspace(           min(Rs_rand(:,3)) , max(Rs_rand(:,3))+1000*10^3,100);
% 
% 
% [x,y]=meshgrid(x,y);
% % syms x y z;
% ct=sqrt((Rs_rand(1,1)-x).^2+(Rs_rand(1,2)-y).^2+(Rs_rand(1,3)-z).^2);
% eqn2= Rs_rand(2,3)-sqrt( (c*delT(2) + ct).^2 - (Rs_rand(2,1)-x).^2 - (Rs_rand(2,2)-y).^2 );
% eqn3= c*delT(3) + ct - sqrt( (Rs_rand(3,1)-x).^2 + (Rs_rand(3,2)-y).^2 + (Rs_rand(3,3)-z).^2);
% eqn4= c*delT(4) + ct - sqrt( (Rs_rand(4,1)-x).^2 + (Rs_rand(4,2)-y).^2 + (Rs_rand(4,3)-z).^2);
% figure
% % earth_sphere(50,'m')
% hold on
% surf(x,y,eqn2)
% % surf(x,y,eqn3)
% % surf(x,y,eqn4)
% % [Solx,Soly,Solz]=solve(eqn2,eqn3,eqn4,x,y,z);
% if isnan(Solx)
%     X=0;
%     Y=0;
%     Z=0;
% 
% elseif length(Solx)>1
%     X1=double(Solx(1));
%     Y1=double(Soly(1));
%     Z1=double(Solz(1));
%     X2=double(Solx(2));
%     Y2=double(Soly(2));
%     Z2=double(Solz(2));
%     if norm([X1,Y1,Z1])>6374*10^(3)&& imag(X1+Y1+Z1)==0
%         X=X1;
%         Y=Y1;
%         Z=Z1;
%     else
%         X=X2;
%         Y=Y2;
%         Z=Z2;
%     end
% else
%    X=double(Solx(1));
%    Y=double(Soly(1));
%    Z=double(Solz(1)); 
% end
% % [X,Y,Z]-R
% if imag(X+Y+Z)~=0
%     sig=sign(real([X,Y,Z]));
%    tmp= sig.*sqrt(conj([X,Y,Z]).*[X,Y,Z]);
%    X=tmp(1);
%    Y=tmp(2);
%    Z=tmp(3);
% end
% r=(norm([X,Y,Z])-6378.14*10^(3))*10^(-3);
% err_r=norm([X,Y,Z]-R)*10^(-3);
% Err_X=((X-6378.14*10^(3))-(R(1)-6378.14*10^(3)))/(R(1)-6378.14*10^(3))*100;
% Err_Y=((Y-6378.14*10^(3))-(R(2)-6378.14*10^(3)))/(R(2)-6378.14*10^(3))*100;
% Err_Z=((Z-6378.14*10^(3))-(R(3)-6378.14*10^(3)))/(R(3)-6378.14*10^(3))*100;
% err=[r,err_r,Err_X,Err_Y,Err_Z];
% end
function [xx,yy,zz] = earth_sphere(varargin)
%EARTH_SPHERE Generate an earth-sized sphere.
%   [X,Y,Z] = EARTH_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURFACE(X,Y,Z) produces a sphere equal to 
%   the radius of the earth in kilometers. The continents will be
%   displayed.
%
%   [X,Y,Z] = EARTH_SPHERE uses N = 50.
%
%   EARTH_SPHERE(N) and just EARTH_SPHERE graph the earth as a 
%   SURFACE and do not return anything.
%
%   EARTH_SPHERE(N,'mile') graphs the earth with miles as the unit rather
%   than kilometers. Other valid inputs are 'ft' 'm' 'nm' 'miles' and 'AU'
%   for feet, meters, nautical miles, miles, and astronomical units
%   respectively.
%
%   EARTH_SPHERE(AX,...) plots into AX instead of GCA.
% 
%  Examples: 
%    earth_sphere('nm') produces an earth-sized sphere in nautical miles
%
%    earth_sphere(10,'AU') produces 10 point mesh of the Earth in
%    astronomical units
%
%    h1 = gca;
%    earth_sphere(h1,'mile')
%    hold on
%    plot3(x,y,z)
%      produces the Earth in miles on axis h1 and plots a trajectory from
%      variables x, y, and z
%   Clay M. Thompson 4-24-1991, CBM 8-21-92.
%   Will Campbell, 3-30-2010
%   Copyright 1984-2010 The MathWorks, Inc. 
%% Input Handling
[cax,args,nargs] = axescheck(varargin{:}); % Parse possible Axes input
error(nargchk(0,2,nargs)); % Ensure there are a valid number of inputs
% Handle remaining inputs.
% Should have 0 or 1 string input, 0 or 1 numeric input
j = 0;
k = 0;
n = 50; % default value
units = 'km'; % default value
for i = 1:nargs
    if ischar(args{i})
        units = args{i};
        j = j+1;
    elseif isnumeric(args{i})
        n = args{i};
        k = k+1;
    end
end
if j > 1 || k > 1
    error('Invalid input types')
end
%% Calculations
% Scale factors
Scale = {'km' 'm'  'mile'            'miles'           'nm'              'au'                 'ft';
         1    1000 0.621371192237334 0.621371192237334 0.539956803455724 6.6845871226706e-009 3280.839895};
% Identify which scale to use
try
    myscale = 6378.1363*Scale{2,strcmpi(Scale(1,:),units)};
catch %#ok<*CTCH>
    error('Invalid units requested. Please use m, km, ft, mile, miles, nm, or AU')
end
     
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = myscale*cosphi*cos(theta);
y = myscale*cosphi*sintheta;
z = myscale*sin(phi)*ones(1,n+1);
%% Plotting
if nargout == 0
    cax = newplot(cax);
    % Load and define topographic data
    load('topo.mat','topo','topomap1');
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    topo2 = [topo(:,181:360) topo(:,1:180)]; %#ok<NODEF>
    
    % Define surface settings
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo2;
    % Create the sphere with Earth topography and adjust colormap
    surface(x,y,z,props,'parent',cax)
    colormap(topomap1)
% Replace the calls to surface and colormap with these lines if you do 
% not want the Earth's topography displayed.
%     surf(x,y,z,'parent',cax)
%     shading flat
%     colormap gray
    
    % Refine figure
    axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])
    view(127.5,30)
else
    xx = x; yy = y; zz = z;
end
end