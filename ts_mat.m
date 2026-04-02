
function [] = ts_mat(expt,savefile)

close all

S = load(['D:\Matlab\' expt '.mat'],'-mat');  %%% loads in data from .mat file into structure S
names = fieldnames(S); %%% gives names of all fields in S. Are 2 for each sweep, one for voltage and one for current; each contains time array and data array

V = getfield(S,names{2});  %% reads in voltage data for this sweep (time & volts)
A = getfield(S,names{1}); %% reads in current data for this sweep (time & amps)
t = V(:,1); %%% time#1 array
v = V(:,2); %%% voltage array
t2 = A(:,1); %%% time#2 array
i = A(:,2); %%% current array
  

tstep = mean(diff(t))*1000000;

plot(t,i,'k-')
axis square

%%%%%%%%%%%% series resistance

pk1j = find(t<0.02);
pk1_t = t(pk1j);
pk1_i = i(pk1j);
[pk1,pk1k] = min(pk1_i);
pk1t = pk1_t(pk1k);

hold on
plot(pk1t,pk1,'m*')

basej = find((t>0.008)&(t<0.01));
base_t = t(basej);
base_i = i(basej);
base = mean(base_i);
disp(base*1000000000000)
plot([base_t(1) base_t(length(base_t))],[base base],'-m');

amp1 = abs(pk1 - base);

vbase = mean(v(basej));
vstep = abs(v(pk1k) - vbase);

rs = (vstep / amp1) / 1000000

%%%%%%%%% input resistance

platj = find ((t>0.019)&(t<0.0198));
plat_t = t(platj);
plat_i = i(platj);
plat = mean(plat_i);

plot([plat_t(1) plat_t(length(plat_t))],[plat plat],'-m');

amp2 = abs(plat - base);

rm = (vstep / amp2) / 1000000

title([' Rs = ' num2str(rs) ' Rm = ' num2str(rm)])
axis([0.005 0.02 pk1-(0.1*amp1) base+(0.1*amp1)])


%%%%%%% membrane capacitance1 - area measure

areabox = abs(0.01 * amp2);

areaj = find((t>0.01)&(t<0.02));
areat = t(areaj);
areai = i(areaj)-base;
area = 0;
for n = 1:((length(areai)-1))
    area = area + (areai(n).*(areat(n+1)-areat(n))) + (0.5*((areai(n+1)-areai(n))*(areat(n+1)-areat(n))));
end

area = abs(area) - areabox;

cm = (area ./ vstep) * 1000000000000



%%%%%%% writing data to file

if nargin > 1 & savefile
    fid = fopen(savefile,'at');
    fprintf(fid,'%s\t %6.3f\t %6.3f\t %6.3f\t %6.3f\n', expt, base*1000000000000, rs, rm, cm);
    fclose(fid);
end

%close all



