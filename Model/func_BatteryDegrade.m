function cap = func_BatteryDegrade(t_hour)

time = t_hour/24;

E_P = 0.38;
c_rate=1/E_P;

a = 8.61e-6;
b = -5.13e-3;
c = 7.63e-1;
d = -6.7e-3;
e = 2.35;
f = 14876;
Ir = c_rate;
Ea = 24.5*10^3;
R = 8.314;
T = 25+273;

load aug15jul16.mat

r_eff=0.9;
Rinc=0.04/12; % resistance increase per month
EOL_eff=1-(1-r_eff)*(Rinc+1);
eff_loss=(r_eff-EOL_eff)/24/30 %monthly efficiency loss divided so that it is once every hour
r_eff=r_eff-(eff_loss*t_hour)
ch_eff=sqrt(r_eff)
dis_eff=sqrt(r_eff)

t_year = length(ACE)-1;
tplot = (0:2:t_year*2)'/3600;

ah=0;

top=sum(max(ACE,0))*(1/dis_eff);
bot=-sum(min(ACE,0))*ch_eff;

    if top>bot
        Cyclesperyear = trapz(tplot,max(ACE,0)/0.38*(1/dis_eff));
    else
        Cyclesperyear = -trapz(tplot,min(ACE,0)/0.38*ch_eff);
    end

    ah = -c_rate*Cyclesperyear*time/(12*30);
    
    qloss = (a*T^2 + b*T + c)*exp((d*T+e)*Ir).*ah + f*time.^(1/2)*exp(-Ea/(R*T));
        
    
cap = 1-qloss/100;