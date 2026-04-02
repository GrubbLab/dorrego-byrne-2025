%%% takes in abf file ('expt') of 10ms current injection steps of increasing
%%% amplitude, plots phase plane plots and determines spike waveform
%%% measures.  Default takes measures from first spike fired at current
%%% threshold; alternatively can choose to take measures from sweep defined
%%% in 'sweep' input.  By defining 'txt' input, can choose to output AP
%%% trace to txt file for analysis in Maxim Volgushev's 'FitAPs' script


function [] = AP_10ms_mat(expt,sweep,txt)

close all

S = load([expt '.mat'],'-mat');  %%% loads in data from .mat file into structure S
names = fieldnames(S); %%% gives names of all fields in S. Are 2 for each sweep, one for voltage and one for current; each contains time array and data array
nsweeps = (length(names)/2); %% so no of sweeps is half no of fields

for j = 1:nsweeps
        
    V = getfield(S,names{(2*j)-1});  %% reads in voltage data for this sweep (time & volts)
    A = getfield(S,names{2*j}); %% reads in current data for this sweep (time & amps)
    t{j} = V(:,1)*1000; %%% time#1 array
    v{j} = V(:,2)*1000; %%% voltage array
    t2{j} = A(:,1)*1000; %%% time#2 array
    i{j} = A(:,2)*1000000000000; %%% current array
    
    %%% then all analysis on these 4 arrays is exactly the same as from
    %%% .asc files
    
    plot(t{j},v{j})
    hold on
    %pause
   
    
    
    ihold(j) = mean(i{j}(1:1950));
    istep(j) = (mean(i{j}(2050:2950)) - ihold(j));
    tstep(j) = mean(diff(t{j}));
    
    vhold(j) = mean(v{j}(1:1950));
    disp(['Trace ' num2str(j) ' I ' num2str(istep(j)) ' Vhold ' num2str(vhold(j))])    
    
    
    [vmax(j) vmax_i(j)] = max(v{j});
    
    if vmax(j)>0
        disp('spike!')
        spikes(j) = 1;
    else
        spikes(j) = 0;
    end
    
    
end

tstep = mean(tstep);
%disp(tstep.*1000)

%%%%%%%%% spike waveform measures from 1st AP

if nargin>1
    j = sweep;
else
    j = find(spikes>0); j = j(1);
end
ts = t{j};
vs = v{j};

figure(2)
subplot(3,1,1)
plot(ts,vs)
hold on

[vmax vmaxi] = max(vs);
vmaxt = ts(vmaxi);
plot(vmaxt,vmax,'g*')

    %%%%%%%%%% finding threshold & other measures
    
window = 20; %%% size of sliding window
a = 1; b = ones(1,window)./window;
vf = filter(b,a,vs); %%% sliding mean
vf = vf(100:length(vf));

    %%%%% saving voltage trace to txt file for FitAPs analysis (Volgushev fits)


    if nargin>2
        savefile = [expt(1:(length(expt)-4)) '_' num2str(j) '.txt'];
        fid = fopen(savefile,'at');
        for n = 1:length(vf)
            fprintf(fid,'%f\n',vf(n));
        end
        fclose(fid);
    end


    %%%%%%%

ts = ts(100:length(ts));
tf = ts-(tstep.*(window-1)/2);
vs = vs(100:length(vs));
plot(tf,vf,'m-')
    
dv = diff(vf);
dvdt = dv ./ tstep;   %%% 1st derivative of smoothed trace
t1 = tf(1:length(tf)-1);    
subplot(3,1,2)
plot(t1,dvdt)
xlabel('s');ylabel('V/s')
hold on    
vf2 = filter(b,a,dvdt);
plot(t1-(tstep.*(window-1)/2),vf2,'m-')
[mx_dvdt mxi_dvdt] = max(dvdt);
tmx_dvdt = t1(mxi_dvdt);    
title(['max dV/dt ' num2str(mx_dvdt)])
plot(tmx_dvdt,mx_dvdt,'*c')
tsh_i = find((t1<tmx_dvdt)&(dvdt<10));    %%% so threshold where 1st derivative goes above 10mV/ms
tshi = tsh_i(length(tsh_i));
plot(t1(tshi),10,'g*')
t_thresh = ts(tshi);    
v_thresh = vs(tshi);
subplot(3,1,1)
plot(t_thresh,v_thresh,'g*')
hh = v_thresh+ (0.5.*(vmax-v_thresh));
hhv1i = find((ts>t_thresh)&(vs>=hh)); hh1i = hhv1i(1);
hhv2i = find((ts>vmaxt)&(vs>=hh));
hh2i = hhv2i(length(hhv2i));
hhv1t = ts(hh1i);
hhv2t = ts(hh2i);
plot([hhv1t hhv2t],[hh hh],'g-')
whh = hhv2t-hhv1t;
title(['vthresh ' num2str(v_thresh) '  vmax ' num2str(vmax) '  whh ' num2str(whh)])
hold off

%%%% phase plot, with optional alternative Vth, and onset rapidness measures

cuti = find(t1>15);  %%% so dispensing with initial transient
v_pp = vf(cuti);
dvdt_pp = dvdt(cuti);

figure(3)
plot(v_pp,dvdt_pp)
axis square
hold on


tsh_i = find(dvdt_pp>=10);
okcount = 0;
for k = 1:(length(tsh_i)-10)
    if min(dvdt_pp(tsh_i(k):tsh_i(k)+10))>=10
        okcount = okcount+1;
        tishi(okcount) = tsh_i(k);
    end
end
tshi = tishi(1);
vth_pp = v_pp(tshi);    %%% is optional vthresh calculated from phase plane plot
dvdt_vth_pp = dvdt_pp(tshi);
%plot(vth_pp,dvdt_vth_pp,'r*')
      
or_v = v_pp(tshi-10:tshi+10);
or_dvdt = dvdt_pp(tshi-10:tshi+10);  %%% 10 points before and after Vthresh point
plot(or_v,or_dvdt,'ro')
p = polyfit(or_v,or_dvdt,1);   %% fitting line to data
yfit = p(1)*([vth_pp-5 vth_pp+5]) + p(2);
plot([vth_pp-5 vth_pp+5],yfit,'r-')
or1 = p(1);
title(['OR = ' num2str(or1)])


%%%% output

first_spike_properties = [j istep(j) vhold(j) v_thresh vmax vmax-v_thresh whh mx_dvdt or1]'
    


%%%%% 2nd derivative measures (optional)



deriv2 = 1;     %% set >0 if 2nd derivative measures desired
if deriv2 >0
 
%%%%% d2V/dt2 iAP

d2V = diff(vf2);
d2Vdt2 = d2V./tstep;
t_d2Vdt2 = t1(1:length(t1)-1);
figure(2)
subplot(3,1,3)
plot(t_d2Vdt2,d2Vdt2,'b-')
title(['d2V/dt2 iAP'])
xlabel('ms');ylabel('V2/s')
hold on

%%%% 2nd derivative measures

%basei = find((t_d2Vdt2 > -1.5)&(t_d2Vdt2 > -.5));
%base = mean(d2Vdt2(basei));
%plot([-1.5 -.5],[base base],'m-')

disp('set zoom for 2nd derivative plot, then press any key')
zoom;pause;zoom;
disp('now click on approximate x locations of two peaks')
[x1,y1] = ginput(1);  %%% plot rough location of 1st peak
pk1_i = find((t_d2Vdt2 > (x1-0.05)) & (t_d2Vdt2 < (x1+0.05)));
d2Vdt2_pk1 = d2Vdt2(pk1_i);
t_d2Vdt2_pk1 = t_d2Vdt2(pk1_i);
[max1 max1_i] = max(d2Vdt2_pk1);
tmax1 = t_d2Vdt2_pk1(max1_i);
plot(tmax1,max1,'mo')

[x2,y2] = ginput(1);  %%% plot rough location of 2nd peak
pk2_i = find((t_d2Vdt2 > (x2-0.05)) & (t_d2Vdt2 < (x2+0.05)));
d2Vdt2_pk2 = d2Vdt2(pk2_i);
t_d2Vdt2_pk2 = t_d2Vdt2(pk2_i);
[max2 max2_i] = max(d2Vdt2_pk2);
tmax2 = t_d2Vdt2_pk2(max2_i);
plot(tmax2,max2,'mo')

amp1 = max1;% - base;
amp1_5 = (0.05 * amp1);% +base  %%% 5% of 1st peak amplitude, a la Kress et al. 08
amp1_5_i = find(((t_d2Vdt2 < (tmax1)) & (d2Vdt2 < (amp1_5))));
amp1_5_i = amp1_5_i(length(amp1_5_i));
t_amp1_5 = t_d2Vdt2(amp1_5_i);
plot(t_amp1_5, amp1_5, 'mo')

RelT5perc = (t_amp1_5 - t_thresh)*1000;
RelTmax1 = (tmax1 - t_thresh)*1000;
RelTmax2 = (tmax2 - t_thresh)*1000;
TmaxDiff = (tmax2 - tmax1)*1000;
MaxRatio = max2 ./ max1;

second_derivative = [RelT5perc, RelTmax1, RelTmax2, TmaxDiff, max1, max2, MaxRatio]'


end
