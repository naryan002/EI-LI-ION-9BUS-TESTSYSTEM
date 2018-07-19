function mpc = case_yl9bus
%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Joe H. Chow's book, p. 70.

%   MATPOWER
%   $Id: case9.m 1559 2010-03-10 18:08:32Z ray $

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
% mpc.bus = [
% 	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
% 	2	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
% 	3	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
% 	4	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
% 	5	1	350	30	0	0	1	1	0	345	1	1.1	0.9;
% 	6	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
% 	7	2	500	35	0	0	1	1	0	345	1	1.1	0.9;
% 	8	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
% 	9	1	150	50	0	0	1	1	0	345	1	1.1	0.9;
% ];

mpc.bus = [
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	3	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	4	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	5	1	143	30	0	0	1	1	0	345	1	1.1	0.9;
	6	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	7	1	159	35	0	0	1	1	0	345	1	1.1	0.9;
	8	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	9	1	198	50	0	0	1	1	0	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    1	200	0	400	    -280  	1    400	     1     	 340   136	0	0	0	0	0	0	0	0	0	0	0;
    2	200	0	400	    -280  	1    400	     1     	 300   120	0	0	0	0	0	0	0	0	0	0	0;
    3	400	0	400	    -280 	1.05 400	     1	     80    32	0	0	0	0	0	0	0	0	0	0	0; 
    6	400	0	400	    -280  	1.05 400	     1     	 80    32	0	0	0	0	0	0	0	0	0	0	0; 
    8	400	0	400	    -280  	1.05 400	     1     	 0      0	0	0	0	0	0	0	0	0	0	0	0; % Wind/solar
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0	    0.0576	0	250	250	250	0	0	1	-360	360;
	4	5	0.017	0.092	0.158	250	250	250	0	0	1	-360	360;
	5	6	0.039	0.17	0.358	150	150	150	0	0	1	-360	360;
	3	6	0	    0.0586	0	300	300	300	0	0	1	-360	360;
	6	7	0.0119	0.1008	0.209	150	150	150	0	0	1	-360	360;
	7	8	0.0085	0.072	0.149	250	250	250	0	0	1	-360	360;
	8	2	0	    0.0625	0	250	250	250	0	0	1	-360	360;
	8	9	0.032	0.161	0.306	250	250	250	0	0	1	-360	360;
	9	4	0.01	0.085	0.176	250	250	250	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% area data
%	area	refbus
mpc.areas = [
	1	5;
];

%% generator cost data
%	1	startup	shutdown	number of cost coeff.	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
	2	3000	0	3	0.1225	1	335;
    2	3000	0	3	0.1 	1	300;
    2	3000	0	3	0       0.1	1;
];
