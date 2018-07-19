clear all

% Runs simulation for the lifetime of the battery

%% Tunables 

life=50400;
Endlife_cap = 0.5;
Starting_cap = 1;
Simu_interval = 720; % one month for degredation interval

ESCap = 20; % MW
ESallowed = 0.8;

k_dollar_coal = 5.00; % dollar/mmbtu, 2014: 2.37, 2017: 2.22
k_dollar_gas = 2.37; % dollar/mmbtu, 2014: 5.00, 2017: 3.67

%emissions factors
k_coal = 206*0.453592; % CO2 kg/mmbtu (0.453592 lbs to kg) 210
k_gas = 117*0.453592; % CO2 kg/mmbtu (0.453592 lbs to kg)
k_coal_nox = 0.227*0.453592;
k_gas_nox = 0.116*0.453592;
k_coal_so2 = 0.15713*0.453592; %0.16163
k_gas_so2 = 0.001*0.453592;

k_co2=[k_gas,k_coal,k_gas,k_gas]; % for heavy coal
k_nox=[k_gas_nox,k_coal_nox,k_gas_nox,k_gas_nox]; % for heavy coal
k_so2=[k_gas_so2,k_coal_so2,k_gas_so2,k_gas_so2]; % for heavy coal

k_dollar=[k_dollar_gas,k_dollar_coal,k_dollar_gas,k_dollar_gas]; % for heavy coal

Flag_wind=0; % 0 for no renewables 1 for renewables
Flag_solar=1;
t_vec = [0:23]';
load_profile = 100*sin(linspace(-0.8*pi,1.2*pi,length(t_vec)))+600;

%wind bellow
wind = Flag_wind*[141.1454209 164.3588944 117.67492 81.14462355 60.24906746 75.97386229 66.24155752 53.00033944 54.0897504 56.28039471 45.30936663 31.58459838 38.24235306 55.3985323 67.19698144 63.89313544 76.75370188 82.61862694 100.3076893 121.2392963 142.839525 148.6643123 161.1639311 145.8123442]; % wind
%solar bellow
solar = Flag_solar*[0 0 0 0 0 0 0 0 2.761714555 187.1547077 307.7808614 260.2846929 282.7520984 260.9285882 373.699438 294.5867518 181.2343721 0 0 0 0 0 0 0];

wind_profile=(wind+solar);

figure
plot(t_vec,load_profile,t_vec,wind_profile)
legend('load','renewable')
fname = 'C:/filepath';
saveas(gca, fullfile(fname,'load'),'fig');

%%
up_min_vec = [4,12,4.5,5];%  heavy coal [12,12,4,4] heavy gas [4,12,4,4]
down_min_vec = [4,12,4.5,5];%  heavy coal [12,12,4,4] heavy gas [4,12,4,4]

alpha_load = 0.03;
alpha_wind = 0.05;

mpc = loadcase('case_yl9bus');
baseMVA = mpc.baseMVA;

% Change line constraints
mpc.branch(:,6) = 1000*ones(size(mpc.branch,1),1);%#of rows,adding one to that and multiplying by 1000 then changing the rateA column to that value

% Change generator limits (ineq, linear)
mpc.gen(1:4,9:10) = [250,100;160,64;240,96;150,60]; % 80% gas, 20% coal another option [160,64;250,100;250,100;140,56] [160,64;275,110;275,110;90,36]

%(rows, columns)

%% Calculating coal cost curves
coef_Hunter_3 = [0.0144,21.077,9080,117483];  % Input/output curve, y(1000 Btu/h), x(MW)
x_Hunter = (10:1:100)';
y_Hunter = polyval(coef_Hunter_3,x_Hunter);

coef_Hunter_2 = polyfit(x_Hunter,y_Hunter,2);
y_Hunter_fit = polyval(coef_Hunter_2,x_Hunter);

coef_Hunter_cost = coef_Hunter_2*k_dollar_coal*10^(-3);
cost_Hunter = polyval(coef_Hunter_cost,x_Hunter);

Hunter_max = 100;

g2_max=mpc.gen(2,9);

g2_ratio=g2_max/Hunter_max;

g2_cost=[coef_Hunter_cost(1)/g2_ratio,coef_Hunter_cost(2),coef_Hunter_cost(3)*g2_ratio];

%% Calculating gas cost curves
coef_Moss_3 = [-0.0013,2.955,6561.2,662025];  % Input/output curve, y(1000 Btu/h), x(MW)
x_Moss = (10:10:750)'; % MW
y_Moss = polyval(coef_Moss_3,x_Moss); % 1000 Btu/h

coef_Moss_2 = polyfit(x_Moss,y_Moss,2); %fits a polynomial line of 2 degrees (Ax^2+Bx+C)
y_Moss_fit = polyval(coef_Moss_2,x_Moss); %uses the two polynominal line to create a new line

coef_Moss_cost = coef_Moss_2*k_dollar_gas*10^(-3);
cost_Moss = polyval(coef_Moss_cost,x_Moss);

Moss_max = 750;

g1_max=mpc.gen(1,9);
g3_max=mpc.gen(3,9);
g4_max=mpc.gen(4,9);

g1_ratio=g1_max/Moss_max;
g3_ratio=g3_max/Moss_max;
g4_ratio=g4_max/Moss_max;

g1_cost = [coef_Moss_cost(1)/g1_ratio,coef_Moss_cost(2),coef_Moss_cost(3)*g1_ratio];
g3_cost = [coef_Moss_cost(1)/g3_ratio,coef_Moss_cost(2),coef_Moss_cost(3)*g3_ratio];
g4_cost = [coef_Moss_cost(1)/g4_ratio,coef_Moss_cost(2),coef_Moss_cost(3)*g4_ratio];

%%
mpc.gencost(1,5:7) = g1_cost; 
mpc.gencost(2,5:7) = g2_cost;
mpc.gencost(3,5:7) = g3_cost;
mpc.gencost(4,5:7) = g4_cost;

%%
% Change reserve cost
Rcost = [5;10;4;4;5;10;4;4]*baseMVA; %heavy coal case

c_on = [3.99*mpc.gen(1,9)+3.67*mpc.gen(1,9)*k_dollar_gas;5.61*mpc.gen(2,9)+7.5*mpc.gen(2,9)*k_dollar_coal;
    3.99*mpc.gen(3,9)+3.67*mpc.gen(3,9)*k_dollar_gas;3.99*mpc.gen(4,9)+3.67*mpc.gen(4,9)*k_dollar_gas]; % $ per occurance

c_off = [3.99*mpc.gen(1,9);5.61*mpc.gen(2,9);3.99*mpc.gen(3,9);3.99*mpc.gen(4,9)]; % $ per occurance
c_rES = 0.1;

mileage = 7.364;
eff=0.9;
Rinc=0.04;
EOL_eff=1-(1-eff)*(Rinc+1);
eff_loss=(eff-EOL_eff)/life*24*30; %13680 is the life of the battery

%% prepare

Current_cap = Starting_cap;
Current_t = 0;
tvec = [];
co2_vec = [];
so2_vec = [];
nox_vec = [];
cap_vec = [];
coal_vec=[];
gas_vec=[];
F_coal=[];
F_gas=[];
g1_nES=[];
g2_nES=[];
g3_nES=[];
g4_nES=[];
g1=[];
g2=[];
g3=[];
g4=[];
g1_vec=[];
g2_vec=[];
g3_vec=[];
g4_vec=[];
BaseLoad = 143+159+198; 

n_interval = length(t_vec);
n_generator = 4;
n_wind = 1;
n_reserve = 4;
n_ES = 1;
n_branch = 9;
n_bus = 9;

%% cost
cost_coef = mpc.gencost(1:4,5:7);
Q_gen = diag(cost_coef(:,1))*baseMVA^2; %making cost in 1st column of cost_coef into the diagonal of a 4x4
c_gen = cost_coef(:,2)*baseMVA;
c_gen_c = cost_coef(:,3);
c_reserve = Rcost(1:n_reserve)*2;

CrES = zeros(n_bus,1);
CrES(8,1) = 1; % the first one is 0.1 and the last is 1 

% power flow (eq, linear)
Bmatrix = zeros(n_bus,n_bus);
branch_info = mpc.branch;
for i = 1:n_branch
    ind_i = branch_info(i,1);
    ind_j = branch_info(i,2);
    bij = 1/branch_info(i,4);
    Bmatrix(ind_i,ind_j) = -bij;
    Bmatrix(ind_j,ind_i) = -bij;
end

Bmatrix = Bmatrix + diag(-sum(Bmatrix,2));

bus_info = mpc.bus;

Cg = zeros(n_bus,n_generator);
Cg(1,1) = 1;
Cg(2,2) = 1;
Cg(3,3) = 1;
Cg(6,4) = 1;

Cw = zeros(n_bus,n_wind);
Cw(8,1) = 1;

% generator limits (ineq, linear)
gen_limits = mpc.gen(1:4,9:10)/baseMVA;

% line constraints (ineq, linear)
Smax = mpc.branch(:,6)/baseMVA;

Bline = zeros(n_branch,n_branch);
for i = 1:n_branch
    ind_i = branch_info(i,1);
    ind_j = branch_info(i,2);
    bij = 1/branch_info(i,4);
    Bline(i,ind_i) = bij;
    Bline(i,ind_j) = -bij;
end

% prepare for cvx
Q_gen_cvx = [];
c_gen_cvx = [];
c_gen_c_cvx = [];
c_reserve_cvx = [];
c_on_cvx = [];
c_off_cvx = [];
c_rES_cvx = [];
Cg_cvx = [];
Cw_cvx = [];
CrES_cvx = [];
Bmatrix_cvx = [];
Sd_cvx = [];
gen_max_cvx = [];
gen_min_cvx = [];
wind_forecast_cvx = [];
Smax_cvx = [];
Bline_cvx = [];


for i = 1:n_interval %n_interval is 24 hours so large loop only runs for 1 day but deg. in code is for 1 month??
    
    % Change loads
    loadamp = load_profile(i)/BaseLoad;
    mpc.bus(:,3:4) = loadamp*bus_info(:,3:4);
    Sd =mpc.bus(:,3)/baseMVA;
    
    % Change generator bounds
    wind_forecast = wind_profile(i);
    mpc.gen(5,9:10) = [wind_forecast,0]; % wind
    
    % Matrices
    Q_gen_cvx = blkdiag(Q_gen_cvx,Q_gen);
    c_gen_cvx = [c_gen_cvx;c_gen];
    c_gen_c_cvx = [c_gen_c_cvx;c_gen_c];
    c_reserve_cvx = [c_reserve_cvx;c_reserve];
    c_on_cvx = [c_on_cvx;c_on];
    c_off_cvx = [c_off_cvx;c_off];
    c_rES_cvx = [c_rES_cvx;c_rES];
    
    Cg_cvx = blkdiag(Cg_cvx,Cg);
    Cw_cvx = blkdiag(Cw_cvx,Cw);
    CrES_cvx = blkdiag(CrES_cvx,CrES);
    Bmatrix_cvx = blkdiag(Bmatrix_cvx,Bmatrix);
    Sd_cvx = [Sd_cvx;Sd];
    
    gen_max_cvx = [gen_max_cvx;gen_limits(:,1)];
    gen_min_cvx = [gen_min_cvx;gen_limits(:,2)];
    wind_forecast_cvx = [wind_forecast_cvx;wind_forecast];
    
    Smax_cvx = [Smax_cvx;Smax];
    Bline_cvx = blkdiag(Bline_cvx,Bline);
    
end

cvx_begin

variables p(n_generator*n_interval) pw(n_wind*n_interval) theta(n_bus*n_interval) r(n_reserve*n_interval)
variable u(n_generator*n_interval) binary
variable v(n_generator*n_interval) binary
variable w(n_generator*n_interval) binary

minimize (quad_form(p,Q_gen_cvx) + c_gen_cvx'*p + c_reserve_cvx'*r + c_on_cvx'*v + c_off_cvx'*w + c_gen_c_cvx'*u)

subject to
-Cg_cvx*p - Cw_cvx*pw + Bmatrix_cvx*theta + Sd_cvx == 0
p + r <= gen_max_cvx.*u
p - r >= gen_min_cvx.*u
0 <= pw <= wind_forecast_cvx/baseMVA

for i = 1:n_interval
    theta_i = theta((i-1)*n_bus+1:i*n_bus);
    u_i = u((i-1)*n_generator+1:i*n_generator);
    -Smax <= Bline*theta((i-1)*n_bus+1:i*n_bus) <= Smax
    
    theta((i-1)*n_bus+1) == 0
    r_i = r((i-1)*n_reserve+1:i*n_reserve);
    r_demand_i = alpha_wind*pw(i)+alpha_load*load_profile(i)/baseMVA;
    sum(r_i) >= r_demand_i
end

r >= 0

u_0_length = max(max(up_min_vec),max(down_min_vec))+1;
u_0 = ones(n_generator*u_0_length,1); % assume all generator on in the previous period
u_extend = [u_0;u];

for k = 1:n_generator
    up_min = up_min_vec(k);
    down_min = down_min_vec(k);
    for i = 2:(n_interval+u_0_length)
        u_i = u_extend((i-1)*n_generator+k);
        u_iminus = u_extend((i-2)*n_generator+k);
        if (n_interval+u_0_length-i) >= down_min
            for j = 1:down_min-1
                u_j = u_extend((i-1+j)*n_generator+k);
                u_iminus - u_i + u_j <= 1
            end
        else
            for j = 1:(n_interval+u_0_length-i)
                u_j = u_extend((i-1+j)*n_generator+k);
                u_iminus - u_i + u_j <= 1
            end
        end
        if (n_interval+u_0_length-i) >= up_min
            for j = 1:up_min-1
                u_j = u_extend((i-1+j)*n_generator+k);
                -u_iminus + u_i - u_j <= 0
            end
        else
            for j = 1:(n_interval+u_0_length-i)
                u_j = u_extend((i-1+j)*n_generator+k);
                -u_iminus + u_i - u_j <= 0
            end
        end
    end
end

for i = 1:n_interval
    u_i = u((i-1)*n_generator+1:i*n_generator);
    if i == 1
        u_iminus = u_0(end-3:end);
    else
        u_iminus = u((i-2)*n_generator+1:(i-1)*n_generator);
    end
    v_i = v((i-1)*n_generator+1:i*n_generator);
    w_i = w((i-1)*n_generator+1:i*n_generator);
    -u_iminus + u_i - v_i <= 0
    u_iminus - u_i - w_i <= 0
end

cvx_end

p_noES = p;
pw_noES = pw;
r_noES = r;
u_noES = u;
v_noES = v;
w_noES = w;

d=1;

while Current_cap > Endlife_cap
    
    Current_t;
    
    if (Current_t) == 0
        Current_cap;
    else
        Current_cap = Starting_cap-(1-f_BatteryDegrade_Xu(Current_t))
    end
    tvec = [tvec;Current_t];
    cap_vec = [cap_vec;Current_cap];
    Current_t = Current_t + Simu_interval %simulation interval needs to be converted from hours to days
    
    ESavailable = ESCap*Current_cap*ESallowed/2; %%% this should just be ESCap
    maxES = ESavailable/baseMVA;
    
    
    %%
    
    clearvars p pw theta r u v w
    
    cvx_begin
    
    variables p(n_generator*n_interval) pw(n_wind*n_interval) theta(n_bus*n_interval) r(n_reserve*n_interval) rES(n_ES*n_interval)
    variable u(n_generator*n_interval) binary
    variable v(n_generator*n_interval) binary
    variable w(n_generator*n_interval) binary
    
    minimize (quad_form(p,Q_gen_cvx) + c_gen_cvx'*p + c_reserve_cvx'*r + c_on_cvx'*v + c_off_cvx'*w + c_rES_cvx'*rES + c_gen_c_cvx'*u)
    
    beta = mileage/24*(1-eff);
    eff=eff-eff_loss;
    
    subject to
    -Cg_cvx*p - Cw_cvx*pw + Bmatrix_cvx*theta + Sd_cvx + CrES_cvx*rES*beta == 0
    p + r <= gen_max_cvx.*u
    p - r >= gen_min_cvx.*u
    0 <= pw <= wind_forecast_cvx/baseMVA
    
    for i = 1:n_interval
        theta_i = theta((i-1)*n_bus+1:i*n_bus);
        u_i = u((i-1)*n_generator+1:i*n_generator);
        -Smax <= Bline*theta((i-1)*n_bus+1:i*n_bus) <= Smax
        
        theta((i-1)*n_bus+1) == 0
        r_i = r((i-1)*n_reserve+1:i*n_reserve);
        rES_i = rES((i-1)*n_ES+1:i*n_ES);
        r_demand_i = alpha_wind*pw(i)+alpha_load*load_profile(i)/baseMVA;
        sum(r_i) + sum(rES_i) >= r_demand_i
    end
    
    r >= 0
    0 <= rES <= maxES
    
    u_0_length = max(max(up_min_vec),max(down_min_vec))+1;
    u_0 = ones(n_generator*u_0_length,1); % assume all generator on in the previous period
    u_extend = [u_0;u];
    for k = 1:n_generator
        up_min = up_min_vec(k);
        down_min = down_min_vec(k);
        for i = 2:(n_interval+u_0_length)
            u_i = u_extend((i-1)*n_generator+k);
            u_iminus = u_extend((i-2)*n_generator+k);
            if (n_interval+u_0_length-i) >= down_min
                for j = 1:down_min-1
                    u_j = u_extend((i-1+j)*n_generator+k);
                    u_iminus - u_i + u_j <= 1
                end
            else
                for j = 1:(n_interval+u_0_length-i)
                    u_j = u_extend((i-1+j)*n_generator+k);
                    u_iminus - u_i + u_j <= 1;
                end
            end
            if (n_interval+u_0_length-i) >= up_min
                for j = 1:up_min-1
                    u_j = u_extend((i-1+j)*n_generator+k);
                    -u_iminus + u_i - u_j <= 0
                end
            else
                for j = 1:(n_interval+u_0_length-i)
                    u_j = u_extend((i-1+j)*n_generator+k);
                    -u_iminus + u_i - u_j <= 0
                end
            end
        end
    end
    
    for i = 1:n_interval
        u_i = u((i-1)*n_generator+1:i*n_generator);
        if i == 1
            u_iminus = u_0(end-3:end);
        else
            u_iminus = u((i-2)*n_generator+1:(i-1)*n_generator);
        end
        v_i = v((i-1)*n_generator+1:i*n_generator);
        w_i = w((i-1)*n_generator+1:i*n_generator);
        -u_iminus + u_i - v_i <= 0
        u_iminus - u_i - w_i <= 0
    end
    
    cvx_end
    
    
    for i = 1:n_interval
        
        p_i = p_noES((i-1)*n_generator+1:i*n_generator);
        p1_noES(i,1) = p_i(1);
        p2_noES(i,1) = p_i(2);
        p3_noES(i,1) = p_i(3);
        p4_noES(i,1) = p_i(4);
        
        u_i = u_noES((i-1)*n_generator+1:i*n_generator);
        u1_noES(i,1) = u_i(1);
        u2_noES(i,1) = u_i(2);
        u3_noES(i,1) = u_i(3);
        u4_noES(i,1) = u_i(4);
        
        r_i = r_noES((i-1)*n_generator+1:i*n_generator);
        r1_noES(i,1) = r_i(1);
        r2_noES(i,1) = r_i(2);
        r3_noES(i,1) = r_i(3);
        r4_noES(i,1) = r_i(4);
        
    end
    
    for i = 1:n_interval
        
        p_i = p((i-1)*n_generator+1:i*n_generator);
        p1(i,1) = p_i(1);
        p2(i,1) = p_i(2);
        p3(i,1) = p_i(3);
        p4(i,1) = p_i(4);
        
        u_i = u((i-1)*n_generator+1:i*n_generator);
        u1(i,1) = u_i(1);
        u2(i,1) = u_i(2);
        u3(i,1) = u_i(3);
        u4(i,1) = u_i(4);
        
        r_i = r((i-1)*n_generator+1:i*n_generator);
        r1(i,1) = r_i(1);
        r2(i,1) = r_i(2);
        r3(i,1) = r_i(3);
        r4(i,1) = r_i(4);
        
    end
    
    if (Current_t)==720
    figure 
    subplot(3,2,1)
    plot(t_vec,p1_noES*baseMVA,'-o',t_vec,p1*baseMVA,'--x','MarkerSize',3)
    hline = refline([0 100]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 250]);
    hline.Color='k';
    hline.LineStyle='--';
    legend('w/o ES','w ES','Location','northoutside','Orientation','Horizontal')
    xlim([0,23])
    
    ylabel('G1 MW (NG)')
    subplot(3,2,2)
    plot(t_vec,p2_noES*baseMVA,'-o',t_vec,p2*baseMVA,'--x','MarkerSize',3)
    hline = refline([0 64]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 160]);
    hline.Color='k';
    hline.LineStyle='--';
    xlim([0,23])

    ylabel('G2 MW (Coal)')
    subplot(3,2,3)
    plot(t_vec,p3_noES*baseMVA,'-o',t_vec,p3*baseMVA,'--x','MarkerSize',3)
    hline = refline([0 96]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 240]);
    hline.Color='k';
    hline.LineStyle='--';
    xlim([0,23])

    ylabel('G3 MW (NG)')
    subplot(3,2,4)
    plot(t_vec,p4_noES*baseMVA,'-o',t_vec,p4*baseMVA,'--x','MarkerSize',3)
    hline = refline([0 60]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 150]);
    hline.Color='k';
    hline.LineStyle='--';
    xlim([0,23])
    ylabel('G4 MW (NG)')
    
    subplot(3,2,5)
    plot(t_vec,pw_noES*baseMVA,'-o',t_vec,pw*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('Renewable MW')
    xlabel('Hours')
    saveas(gca, fullfile(fname,'Dispatch_Init'),'fig');

    figure %7
    subplot(3,2,1)
    plot(t_vec,r1_noES*baseMVA,'-o',t_vec,r1*baseMVA,'--x','MarkerSize',3)
    legend('w/o ES','w ES','Location','northoutside','Orientation','Horizontal')
    title('Reserve dispatch (MW)')
    xlim([0,23])
%     ylim([0,30])
    ylabel('G1 MW (NG)')
    subplot(3,2,2)
    plot(t_vec,r2_noES*baseMVA,'-o',t_vec,r2*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
%     ylim([0,30])
    ylabel('G2 MW (Coal)')
    subplot(3,2,3)
    plot(t_vec,r3_noES*baseMVA,'-o',t_vec,r3*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
%     ylim([0,30])
    ylabel('G3 MW (NG)')
    subplot(3,2,4)
    plot(t_vec,r4_noES*baseMVA,'-o',t_vec,r4*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
%     ylim([0,30])
    ylabel('G4 MW (NG)')
    subplot(3,2,5)
    plot(t_vec,zeros(length(t_vec),1),'-o',t_vec,rES*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
%     ylim([0,30])
    ylabel('Battery MW')
    xlabel('Hours')
    saveas(gca, fullfile(fname,'Dispatch_Init'),'fig');
    end
    
    co2_total_noES = [];
    so2_total_noES = [];
    nox_total_noES = [];
    cost_p_noES = [];
    cost_r_noES = [];
    cost_total_noES = [];
    gen1=[];
    gen2=[];
    gen3=[];
    gen4=[];
    gen1_noES=[];
    gen2_noES=[];
    gen3_noES=[];
    gen4_noES=[];
    Fuel_coal_noES = [];
    Fuel_gas_noES = [];
    
    for i = 1:length(t_vec)
        
        p_1 = p1_noES(i)*baseMVA;
        p_2 = p2_noES(i)*baseMVA;
        p_3 = p3_noES(i)*baseMVA;
        p_4 = p4_noES(i)*baseMVA;
        
        gen1_nES(i)=p_1;
        gen2_nES(i)=p_2;
        gen3_nES(i)=p_3;
        gen4_nES(i)=p_4;
                
        cost_1 = polyval(mpc.gencost(1,5:7),p_1);
        cost_2 = polyval(mpc.gencost(2,5:7),p_2);
        cost_3 = polyval(mpc.gencost(3,5:7),p_3);
        cost_4 = polyval(mpc.gencost(4,5:7),p_4);
        
        cost_p_noES(i,1) = cost_1 + cost_2 + cost_3 + cost_4;
        
        r_i = r_noES((i-1)*n_reserve+1:i*n_reserve);
        cost_r_noES(i,1) = c_reserve'*r_i;
        
        cost_total_noES(i,1) = cost_p_noES(i,1) + cost_r_noES(i,1);
        
        Fuel1 = cost_1/k_dollar(1,1); %MMBtu
        Fuel2 = cost_2/k_dollar(1,2);
        Fuel3 = cost_3/k_dollar(1,3);
        Fuel4 = cost_4/k_dollar(1,4);
        
        co2_1 = cost_1/k_dollar(1,1)*k_co2(1,1);
        co2_2 = cost_2/k_dollar(1,2)*k_co2(1,2);
        co2_3 = cost_3/k_dollar(1,3)*k_co2(1,3);
        co2_4 = cost_4/k_dollar(1,4)*k_co2(1,4);
        
        so2_1 = cost_1/k_dollar(1,1)*k_so2(1,1);
        so2_2 = cost_2/k_dollar(1,2)*k_so2(1,2);
        so2_3 = cost_3/k_dollar(1,3)*k_so2(1,3);
        so2_4 = cost_4/k_dollar(1,4)*k_so2(1,4);
        
        nox_1 = cost_1/k_dollar(1,1)*k_nox(1,1);
        nox_2 = cost_2/k_dollar(1,2)*k_nox(1,2);
        nox_3 = cost_3/k_dollar(1,3)*k_nox(1,3);
        nox_4 = cost_4/k_dollar(1,4)*k_nox(1,4);
        
        if (p_1<1)
            Fuel1=0;
            co2_1=0;
            so2_1=0;
            nox_1=0;
        end
        if (p_2<1)
            Fuel2=0;
            co2_2=0;
            so2_2=0;
            nox_2=0;
        end
        if (p_3<1)
            Fuel3=0;
            co2_3=0;
            so2_3=0;
            nox_3=0;
        end
        if (p_4<1)
            Fuel4=0;
            co2_4=0;
            so2_4=0;
            nox_4=0;
        end
        
        co2_p = co2_1+co2_2+co2_3+co2_4;
        so2_p = so2_1+so2_2+so2_3+so2_4;
        nox_p = nox_1+nox_2+nox_3+nox_4;
        
        coal_p = Fuel2; % heavy coal
        gas_p = Fuel3 + Fuel4+Fuel1; % heavy coal
        
        co2_total_noES(i,1) = co2_p; %represents a day's worth of emissions
        so2_total_noES(i,1) = so2_p; %represents a day's worth of emissions
        nox_total_noES(i,1) = nox_p; %represents a day's worth of emissions
        Fuel_coal_noES(i,1) = coal_p;
        Fuel_gas_noES(i,1) = gas_p;
        Fuel_G1_noES(i,1) = Fuel1;
        Fuel_G2_noES(i,1) = Fuel2;
        Fuel_G3_noES(i,1) = Fuel3;
        Fuel_G4_noES(i,1) = Fuel4;
    end
    
    g1_nES(d)=sum(gen1_nES);
    g2_nES(d)=sum(gen2_nES);
    g3_nES(d)=sum(gen3_nES);
    g4_nES(d)=sum(gen4_nES);
    f1_nES(d)=sum(Fuel_G1_noES);
    f2_nES(d)=sum(Fuel_G2_noES);
    f3_nES(d)=sum(Fuel_G3_noES);
    f4_nES(d)=sum(Fuel_G4_noES);

    
    co2_total = [];
    so2_total = [];
    nox_total = [];
    cost_p = [];
    cost_r = [];
    cost_total = [];
    Fuel_coal = [];
    Fuel_gas = [];

    
    for i = 1:length(t_vec)
        
        p_1 = p1(i)*baseMVA;
        p_2 = p2(i)*baseMVA;
        p_3 = p3(i)*baseMVA;
        p_4 = p4(i)*baseMVA;
        
        gen1(i)=p_1;
        gen2(i)=p_2;
        gen3(i)=p_3;
        gen4(i)=p_4;
        
        cost_1 = polyval(mpc.gencost(1,5:7),p_1);
        cost_2 = polyval(mpc.gencost(2,5:7),p_2);
        cost_3 = polyval(mpc.gencost(3,5:7),p_3);
        cost_4 = polyval(mpc.gencost(4,5:7),p_4);
        
        cost_p(i,1) = cost_1 + cost_2 + cost_3 + cost_4;
        
        r_i = r((i-1)*n_reserve+1:i*n_reserve);
        cost_r(i,1) = c_reserve'*r_i;
        
        cost_total(i,1) = cost_p(i,1) + cost_r(i,1);
        
        co2_1 = cost_1/k_dollar(1,1)*k_co2(1,1);
        co2_2 = cost_2/k_dollar(1,2)*k_co2(1,2);
        co2_3 = cost_3/k_dollar(1,3)*k_co2(1,3);
        co2_4 = cost_4/k_dollar(1,4)*k_co2(1,4);
        
        so2_1 = cost_1/k_dollar(1,1)*k_so2(1,1);
        so2_2 = cost_2/k_dollar(1,2)*k_so2(1,2);
        so2_3 = cost_3/k_dollar(1,3)*k_so2(1,3);
        so2_4 = cost_4/k_dollar(1,4)*k_so2(1,4);
        
        nox_1 = cost_1/k_dollar(1,1)*k_nox(1,1);
        nox_2 = cost_2/k_dollar(1,2)*k_nox(1,2);
        nox_3 = cost_3/k_dollar(1,3)*k_nox(1,3);
        nox_4 = cost_4/k_dollar(1,4)*k_nox(1,4);
        
        Fuel1 = cost_1/k_dollar(1,1); %MMBtu
        Fuel2 = cost_2/k_dollar(1,2);
        Fuel3 = cost_3/k_dollar(1,3);
        Fuel4 = cost_4/k_dollar(1,4);
        
        if (p_1<1)
            Fuel1=0;
            co2_1=0;
            so2_1=0;
            nox_1=0;
        end
        if (p_2<1)
            Fuel2=0;
            co2_2=0;
            so2_2=0;
            nox_2=0;
        end
        if (p_3<1)
            Fuel3=0;
            co2_3=0;
            so2_3=0;
            nox_3=0;
        end
        if (p_4<1)
            Fuel4=0;
            co2_4=0;
            so2_4=0;
            nox_4=0;
        end
        
        co2_p = co2_1+co2_2+co2_3+co2_4;
        so2_p = so2_1+so2_2+so2_3+so2_4;
        nox_p = nox_1+nox_2+nox_3+nox_4;
        
        coal_p = Fuel2; % heavy coal
        gas_p = Fuel3 + Fuel4 + Fuel1; % heavy coal

        co2_total(i,1) = co2_p;
        so2_total(i,1) = so2_p;
        nox_total(i,1) = nox_p;
        Fuel_coal(i,1) = coal_p;
        Fuel_gas(i,1) = gas_p;
        Fuel_G1(i,1) = Fuel1;
        Fuel_G2(i,1) = Fuel2;
        Fuel_G3(i,1) = Fuel3;
        Fuel_G4(i,1) = Fuel4;
        
        
    end

    g1(d)=sum(gen1);
    g2(d)=sum(gen2);
    g3(d)=sum(gen3);
    g4(d)=sum(gen4);
    f1(d)=sum(Fuel_G1);
    f2(d)=sum(Fuel_G2);
    f3(d)=sum(Fuel_G3);
    f4(d)=sum(Fuel_G4);
    d=d+1;
   
    F_coal_delta = sum(Fuel_coal) - sum(Fuel_coal_noES); %total for one day in each month
    F_gas_delta = sum(Fuel_gas) - sum(Fuel_gas_noES); %total for one day in each month
    co2_delta = sum(co2_total) - sum(co2_total_noES);
    so2_delta = sum(so2_total) - sum(so2_total_noES);
    nox_delta = sum(nox_total) - sum(nox_total_noES);
    g1_delta=sum(g1)-sum(g1_nES);
    g2_delta=sum(g2)-sum(g2_nES);
    g3_delta=sum(g3)-sum(g3_nES);
    g4_delta = sum(g4)-sum(g4_nES);
    
    coal_vec = [coal_vec;F_coal_delta];
    gas_vec = [gas_vec;F_gas_delta];
    co2_vec = [co2_vec;co2_delta];
    so2_vec = [so2_vec;so2_delta];
    nox_vec = [nox_vec;nox_delta];% possitive if ES increased emissions
    g1_vec = [g1_vec;g1_delta];
    g2_vec = [g2_vec;g2_delta];
    g3_vec = [g3_vec;g3_delta];
    g4_vec = [g4_vec;g4_delta];

    curt_ES = sum(pw*baseMVA)/sum(wind_profile)
    curt_noES = sum(pw_noES*baseMVA)/sum(wind_profile)

if (Current_t)==720
    figure 
    subplot(2,2,1)
    plot(t_vec,Fuel_G1_noES,'-o',t_vec,Fuel_G1,'--x','MarkerSize',3)
    legend('w/o ES','w ES','Location','northoutside','Orientation','Horizontal')
    title('Fuel (MMBtu)')
    xlim([0,23])
    ylabel('G1 MMBtu (Coal)')
    
    subplot(2,2,2)
    plot(t_vec,Fuel_G2_noES,'-o',t_vec,Fuel_G2,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('G2 MMBtu (Coal)')
    
    subplot(2,2,3)
    plot(t_vec,Fuel_G3_noES,'-o',t_vec,Fuel_G3,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('G3 MMBtu (NG)')
    
    subplot(2,2,4)
    plot(t_vec,Fuel_G4_noES,'-o',t_vec,Fuel_G4,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('G4 MMBtu (NG)')
    saveas(gca, fullfile(fname,'Fuel_Init'),'fig');
    end
end
%%

    figure 
    subplot(3,2,1)
    plot(t_vec,p1_noES*baseMVA,'-o',t_vec,p1*baseMVA,'--x','MarkerSize',3)
    legend('w/o ES','w ES','Location','northoutside','Orientation','Horizontal')
    title('Generator dispatch (MW)')
    hline = refline([0 100]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 250]);
    hline.Color='k';
    hline.LineStyle='--';
    xlim([0,23])
    ylabel('G1 MW (NG)')
    
    subplot(3,2,2)
    plot(t_vec,p2_noES*baseMVA,'-o',t_vec,p2*baseMVA,'--x','MarkerSize',3)
    hline = refline([0 64]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 160]);
    hline.Color='k';
    hline.LineStyle='--';
    xlim([0,23])
    ylabel('G2 MW (Coal)')
    
    subplot(3,2,3)
    plot(t_vec,p3_noES*baseMVA,'-o',t_vec,p3*baseMVA,'--x','MarkerSize',3)
    hline = refline([0 96]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 240]);
    hline.Color='k';
    hline.LineStyle='--';
    xlim([0,23])
    ylabel('G3 MW (NG)')
    xlabel('Hours')
    
    subplot(3,2,4)
    plot(t_vec,p4_noES*baseMVA,'-o',t_vec,p4*baseMVA,'--x','MarkerSize',3)
    hline = refline([0 60]);
    hline.Color='k';
    hline.LineStyle='--';
    hline = refline([0 150]);
    hline.Color='k';
    hline.LineStyle='--';
    xlim([0,23])
    ylabel('G4 MW (NG)')
    xlabel('Hours')
    
    subplot(3,2,5)
    plot(t_vec,pw_noES*baseMVA,'-o',t_vec,pw*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('Renewable MW')
    xlabel('Hours')
    saveas(gca, fullfile(fname,'Dispatch_Final'),'fig');
    
    
    figure %7
    subplot(3,2,1)
    plot(t_vec,r1_noES*baseMVA,'-o',t_vec,r1*baseMVA,'--x','MarkerSize',3)
    legend('w/o ES','w ES','Location','northoutside','Orientation','Horizontal')
    title('Reserve dispatch (MW)')
    xlim([0,23])
    ylabel('G1 MW (NG)')
    
    subplot(3,2,2)
    plot(t_vec,r2_noES*baseMVA,'-o',t_vec,r2*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('G2 MW (Coal)')
    
    subplot(3,2,3)
    plot(t_vec,r3_noES*baseMVA,'-o',t_vec,r3*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('G3 MW (NG)')
    
    subplot(3,2,4)
    plot(t_vec,r4_noES*baseMVA,'-o',t_vec,r4*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('G4 MW (NG)')
    
    subplot(3,2,5)
    plot(t_vec,zeros(length(t_vec),1),'-o',t_vec,rES*baseMVA,'--x','MarkerSize',3)
    xlim([0,23])
    ylabel('Battery MW')
    xlabel('Hours')
    saveas(gca, fullfile(fname,'Reserves_Final'),'fig');
    
total_co2 = sum(30*co2_vec)%Simu_interval is 30*24
total_so2 = sum(30*so2_vec)
total_nox = sum(30*nox_vec)
total_coal = sum(30*coal_vec)
total_gas = sum(30*gas_vec)
total_coal_MW=30*(sum(g2_vec))
total_gas_MW=30*(sum(g3_vec)+sum(g4_vec)+sum(g1_vec))

label={'PriceSwap_NG_S'};
impacts(1,:)=[total_co2,total_so2,total_nox,total_gas,total_coal,total_gas_MW,total_coal_MW];
xlswrite('impacts_Xu.xlsx',label,1,'A13')
xlswrite('impacts_Xu.xlsx',impacts,1,'B13')

%%
figure
plot(tvec,co2_vec,'-o') % tvec is in hours not days
xlabel('Hours')
ylabel('CO_2 (kg)')
saveas(gca, fullfile(fname,'CO2'),'fig');

figure
plot(tvec,coal_vec,'-o',tvec,gas_vec,'-o') % tvec is in hours not days
xlabel('Hours')
ylabel('Fuel (mmBtu)')
legend('Coal','Gas')
saveas(gca, fullfile(fname,'Fuel'),'fig');

figure
plot(tvec/24/30,g1,'-o',tvec/24/30,g2,'-o',tvec/24/30,g3,'-o',tvec/24/30,g4,'-o')
legend('g1','g2','g3','g4')
xlabel('Months')
ylabel('Generation (MWh)')
saveas(gca, fullfile(fname,'Gen_ES'),'fig');

figure
plot(tvec/24/30,g1_nES,'-o',tvec/24/30,g2_nES,'-o',tvec/24/30,g3_nES,'-o',tvec/24/30,g4_nES,'-o')
legend('g1_nES','g2_nES','g3_nES','g4_nES')
xlabel('Months')
ylabel('Generation (MWh)')
saveas(gca, fullfile(fname,'Gen_noES'),'fig');