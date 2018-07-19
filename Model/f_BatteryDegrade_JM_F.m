% JM degredation model seperate from code
function cap = f_BatteryDegrade_JM_F(t_hour)

load aug15jul16.mat

time = t_hour/24;

int=3600/2; %number of intervals in one day
t=length(ACE)-int+1;

r_eff=0.9;
Rinc=0.04/12; % resistance increase per month
EOL_eff=1-(1-r_eff)*(Rinc+1);
eff_loss=(r_eff-EOL_eff)/24/30 %monthly efficiency loss divided so that it is once every hour
r_eff=r_eff-(eff_loss*t_hour)
ch_eff=sqrt(r_eff)
dis_eff=sqrt(r_eff)

a_JM = 0.37;
b_JM = 0.42;
c_JM = 0;
d_JM = 6.5e-3;
e_JM = 6e-3;

f = 14876;
Ea = 24.5*10^3;
R = 8.314;
T = 25+273;

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

sig=1;
fin=length(SOCd);

for h=1:fin %loop for SOCd because it is constant across an hour
    for i = 1:int %loop for each hour because ACE etc. changes
        SOC_i = SOC_m(sig);
        p_discharge = max(ACE(sig),0)/0.38*dis_eff+SOCd(h)*int;
        p_charge = -min(ACE(sig),0)/0.38*(1/ch_eff)+SOCd(h)*int;
        CapLoss(sig,1) = (b_JM*(SOC_i-a_JM)^2+c_JM*p_discharge+d_JM*p_charge+e_JM*p_charge^2)/400/int; 
        sig=sig+1;
    end
end

cap = 1-((sum(CapLoss)*time/30/12)); % t is in days
