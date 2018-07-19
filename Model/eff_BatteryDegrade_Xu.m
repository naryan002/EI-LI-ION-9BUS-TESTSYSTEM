% JM degredation model seperate from code
function cap = f_BatteryDegrade_Xu(t_hour)

load aug15jul16.mat

time = t_hour/24;

int=3600/2; %number of intervals in one day
t=length(ACE)-int+1;

r_eff=1;
%Rinc=0.04/12; % resistance increase per month
%EOL_eff=1-(1-r_eff)*(Rinc+1);
eff_loss=0 %monthly efficiency loss divided so that it is once every hour
%r_eff=r_eff-(eff_loss*t_hour)
ch_eff=sqrt(r_eff)
dis_eff=sqrt(r_eff)

f = 14876;
Ea = 24.5*10^3;
R = 8.314;
T = 25+273;

%Xu model constants
    %general
a_sei=0.0575;
b_sei=121;
    %DOD stress model
k_1=1.4e5;
k_2=-0.501;
k_3=-1.23e5;
    %SOC stress model
SOC_ref=0.5;
k_sig=1.04;
    %time aging
k_t=4.14e-10;

SOC = [];
SOCd=[];
SOC_m=[];
dSOC_m=[];
SOC_c = 0.5;
SOC(1,1)=SOC_c;
hour=1;
E_P = 0.38;
c_rate=1/E_P;

for j=1:int:t  
    for i=1:int % first interval is to determine the needed manipulation of SOC
        p=j+i-1;
        if ACE(p)>0 %discharging
            dSOC(p,1) = -ACE(p)/int/E_P*(1/dis_eff); %this signal is every two seconds
        else %charging
            dSOC(p,1) = -ACE(p)/int/E_P*ch_eff; %this signal is every two seconds
        end
        SOC(p,1) = SOC_c + dSOC(p,1);
        SOC_c = SOC(p,1);
    end
    SOCd(hour)=(0.5-SOC(p,1))/(int);
    SOC_c=SOC(j,1);
    SOC_m(j,1)=SOC(j,1);
    for i=1:int %adjusted the signal accordingly
        p=j+i-1;
        if ACE(p)>0 %discharging
            dSOC_m(p,1) = -(ACE(p)/int/E_P)*(1/dis_eff)+SOCd(hour);
        else %charging
            dSOC_m(p,1) = -(ACE(p)/int/E_P)*(ch_eff)+SOCd(hour);
        end
        SOC_m(p,1) = SOC_c + dSOC_m(p,1);
        if SOC_m(p,1)<0.2
            SOC_m(p,1)=0.2;
        elseif SOC_m(p,1)>0.8
            SOC_m(p,1)=0.8;
        end
        SOC_c = SOC_m(p,1);
    end
    SOC_c=SOC_m(p,1);
    hour=hour+1;
    r_eff=r_eff-eff_loss;
    ch_eff=sqrt(r_eff);
    dis_eff=sqrt(r_eff);
end

% Loop for calculating DOD, average SOC and time degradation for each cycle
j=1;
cy=1;
f_d=0;
SOC_m_avg=0;
DOD=[];
S_DOD=[];
S_SOC=[];
SOC_avg=[];
f_c=[];
rain=1;

t_year = length(ACE)-1;
t_day = (0:30:1500)';
tplot = (0:2:t_year*2)'/3600;
splot = (0:2:length(SOC_m)*2-1);

ext=sig2ext(SOC_m);

rf = rainflow(ext);
DOD_r=rf(1,:)*2;
SOC_r_avg=rf(2,:);
N=rf(3,:);

for i=1:length(DOD_r)
    S_DOD(i)=(k_1*DOD_r(i)^k_2+k_3)^(-1);
    S_SOC(i)=exp(k_sig*(SOC_r_avg(i)-SOC_ref));
    f_c(i)=S_DOD(i)*S_SOC(i)*N(i);
end

%total cycle losses Xu
f_cy=sum(f_c); %(2)
S_SOC_avg=exp(k_sig*(mean(SOC_m)-SOC_ref)); %(25)
S_t=[];
   
S_t=k_t*time*24*3600; %changing number od days to seconds
f_t=S_t*S_SOC_avg; %(21)
f_d=f_t+f_cy*time/360; %(3) f_cy is for the entire year

L = (1-a_sei*exp(-b_sei*f_d)-(1-a_sei)*exp(-f_d)); %Xu

cap = 1-L; % t is in days
