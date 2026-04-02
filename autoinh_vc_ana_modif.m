function [] = autoinh_vc_ana(expt)

disp(expt)
close all

S = load([expt '.mat'],'-mat');  %%% loads in data from .mat file into structure S
names = fieldnames(S); %%% gives names of all fields in S. Are 2 for each sweep, one for voltage and one for current; each contains time array and data array
nsweeps = (length(names)/2); %% so no of sweeps is half no of fields

npulses = 1;   %%% number of pulses
testt = 10;   %%% duration of pulse(s)

start_pulse = 240;  %%% start time of depolarisation
if npulses==2   %%% for paired-pulse experiments
    start_pulse2 = 400;
end

start_test1 = 25;  %% time of start of first test pulse for P4
iti = 35;  %%% inter-test interval (from start of one to start of next)


%%%% generating mean trace

for n=1:nsweeps
    Vee = getfield(S,names{(2*n)-1});  
    Ay = getfield(S,names{2*n}); 
    t(n,:) = Vee(:,1).*1000; %%% time array in ms
    a(n,:) = Vee(:,2).*1000000000000; %%% current array in pA
    v(n,:) = Ay(:,2).*1000; %%% voltage array in mV
    subplot(3,1,1)
    plot(t(n,:),v(n,:))
    subplot(3,1,2)
    plot(t(n,:),a(n,:))
    hold on
    %pause
end

if nsweeps>1
    meana = mean(a);
else
    meana = a(1,:);
end
basei = find(t(1,:)<start_test1);
meana = meana - mean(meana(basei));
subplot(3,1,3)
plot(t(1,:),meana)
hold on

j = find((t(1,:)>=(start_pulse+testt+10))&(t(1,:)<(start_pulse+testt+510)));  %%% Standard analysis, from 10ms after step offset to safely allow for passive changes between expts, and 500ms in duration
t_j = t(1,j);
a_j = meana(j);
plot(t_j,a_j,'b-')

       
a_j_A = a_j./1000000000000;    %%%% current in amps
t_j_s = t_j./1000;     %%% time in seconds
area1_C = -trapz(t_j_s,a_j_A); %%% charge in C
area1_pC = area1_C.*1000000000000  %%%% charge in pC
        


end 



