%%% takes in an image and a file with co-ordinates of pre-detected puncta
%%% measures colocalisation of those puncta with label in a different
%%% channel. Optional 'translate' input is 2x1 array with x+y distances in
%%% um for translating the measurement channel to get estimates of
%%% colocalisation likelihood by chance. Outputs an array of binary values
%%% (cloc), where 1=colocalised and 0=not colocalised for each punctum.
%%% Also an overall colocalised % (perc_cloc) for this image

function [] = syncoloc(image_file,xyz_file,translate)

close all

%%% setting parameters
n_chan = 4;
puncta_chan = 3;    %%% channel with pre-detected puncta
analysis_chan = 2;  %%% channel with signal for colocalisation analysis

xy_res = 0.0329471; %%% XY pixel size in um
umpix = round(1/xy_res); %%% no of pixels per um
z_res = 0.3725992;  %%% Z step size in um

perc = 10; %%% percent of peak to define puncta width
backg_pixn = 20; %%% no of pixels to define background mean (with xy_res = 0.03, 30=1um)
m = 1.5; %%% no of um to view punctum profile over, from punctum centroid in each direction
    
%%% reading in the image file
I = tiffread(image_file);
n_z = length(I)./n_chan; %%% tiffread gives structure I with separate image for each channel+Z position

%%% reading in the puncta file
fid = fopen(xyz_file,'r');
n = 0;
while feof(fid)~=1
    n = n+1;
    x(n) = fscanf(fid,'%f',1);
    y(n) = fscanf(fid,'%f',1);
    z(n) = fscanf(fid,'%f',1);
end
fclose(fid);

%%% working through each pre-detected punctum
for n = 1:length(x) 
    
    disp(['Punctum ' num2str(n)])

    %%%% finding central z-plane for this punctum, and relevant images
    z_int(n) = round(z(n)./z_res) + 1;
    %z_int(n) = 1; %%% switch this line for previous one if image is a MIP/singleZplane
    i(n) = (z_int(n)-1).*n_chan + puncta_chan;  %%% is image number for punctum
    j(n) = (z_int(n)-1).*n_chan + analysis_chan; %%% is image number for analysis

    %%%% plotting punctum position
    subplot(2,3,1)
    Ch1 = I(i(n)).data; %%% specifying the punctum image as 'Ch1'
    mxfi(n) = max(max(Ch1));
    imshow(Ch1,[0 mxfi(n)]) %% punctum full image, scaled to min-max
    title(['Punctum ' num2str(n)])
    colormap('jet')
    hold on
    xpix(n) = x(n)./xy_res;
    ypix(n) = y(n)./xy_res;
    plot(xpix(n),ypix(n),'yo','markerfacecolor','y','markersize',3) %%% overlaying position of current punctum
    hold off

    %%% plotting zoomed in image of punctum in Ch1
    subplot(2,3,2)
    imshow(Ch1,[0 mxfi(n)])
    title('Ch1 zoom')
    colormap('jet')
    axis([xpix(n)-(m*umpix) xpix(n)+(m*umpix) ypix(n)-(m*umpix) ypix(n)+(m*umpix)]) %%% zoomed-in image of punctum +/-2um from centroid
    hold on
    %plot([xpix(n)-umpix xpix(n)+umpix],[ypix(n) ypix(n)],'-r') 
    
    %%% finding fluorescence of line profile through punctum centre
    Ch1_linex{n} = [round(xpix(n)-(m*umpix)):round(xpix(n)+(m*umpix))];
    Ch1_liney{n} = round(ypix(n))*ones(1,length(Ch1_linex{n}));
    for k = 1:length(Ch1_linex{n})
        Ch1_linef{n}(k) = Ch1(round(ypix(n)),Ch1_linex{n}(k)); %%%% matlab uses rows then columns, so this has to be in [y,x] co-ordinates!
        Ch1_linef3{n}(k) = mean([Ch1(round(ypix(n))-1,Ch1_linex{n}(k)) Ch1(round(ypix(n)),Ch1_linex{n}(k)) Ch1(round(ypix(n))+1,Ch1_linex{n}(k))]); %%% mean of 3 pixels centred on punctum y-centre
        Ch1_linef7{n}(k) = mean([Ch1(round(ypix(n))-3,Ch1_linex{n}(k)) Ch1(round(ypix(n))-2,Ch1_linex{n}(k)) Ch1(round(ypix(n))-1,Ch1_linex{n}(k)) Ch1(round(ypix(n)),Ch1_linex{n}(k)) Ch1(round(ypix(n))+1,Ch1_linex{n}(k)) Ch1(round(ypix(n))+2,Ch1_linex{n}(k)) Ch1(round(ypix(n))+3,Ch1_linex{n}(k))]); %%% mean of 7 pixels centred on punctum y-centre
    end
    %plot(Ch1_linex{n},Ch1_liney{n},'y-')
    plot(Ch1_linex{n},Ch1_liney{n}+3,'y-')
    plot(Ch1_linex{n},Ch1_liney{n}-3,'y-')
    subplot(2,3,3)
    %plot(Ch1_linex{n},Ch1_linef{n})
    plot(Ch1_linex{n},Ch1_linef7{n},'g-')
    title('Ch1 profile')
    hold on
    
    %%% finding x-range of punctum using Gaussian fit
    Ch1_backg(n) = mean([Ch1_linef7{n}(1:backg_pixn) Ch1_linef7{n}(end-backg_pixn:end)]);
    plot([Ch1_linex{n}(1) Ch1_linex{n}(backg_pixn)],[Ch1_backg(n) Ch1_backg(n)],'-m')
    plot([Ch1_linex{n}(end-backg_pixn) Ch1_linex{n}(end)],[Ch1_backg(n) Ch1_backg(n)],'-m')
    Ch1_sub{n} = Ch1_linef7{n} - Ch1_backg(n);
    plot(Ch1_linex{n},Ch1_sub{n},'-b')
    xfit = Ch1_linex{n}'; 
    yfit = Ch1_sub{n}';
    fitopt = fitoptions('Method','NonlinearLeastSquares','Lower',[0 (x(n)./xy_res)-20 0],'Upper',[inf (x(n)./xy_res)+20 inf]);  %% setting limits for x position of Gaussian peak (param b) at +-20 pixels from punctum centroid
    fitsub = fit(xfit,yfit,'gauss1',fitopt);
    plot(fitsub)
    legend('off')
    ylabel('F')
    xlabel('pixels')
    params = coeffvalues(fitsub);
    a = params(1);
    b = params(2);
    c = params(3);
    for p = 1:length(xfit)
        fx(p) = a*exp(-((xfit(p)-b)/c)^2);
    end
    %plot(xfit,fx,'mo')
    real_var(n) = var(yfit);
    dif = fx - yfit';
    diff_var(n) = var(dif);
    pe(n) = 100*(1-(diff_var(n)/real_var(n)));
    %title(['pe = ' num2str(pe(n))])
    
    %%%% finding x positions defining puntum width
    thresh = a/100*perc; %%% percent of max for start/end positions
    plot([xfit(1) xfit(end)],[thresh thresh],'-g')
    threshi = find(fx>thresh);
    xt1(n) = xfit(threshi(1)); %% x start of punctum
    xt2(n) = xfit(threshi(end)); %% x end of puncum
    plot([xt1(n) xt1(n)],[min(Ch1_sub{n}) max(Ch1_sub{n})],'g-')
    plot([xt2(n) xt2(n)],[min(Ch1_sub{n}) max(Ch1_sub{n})],'g-')
    hold off

    subplot(2,3,2)
    plot([xt1(n) xt1(n)],[ypix(n)-(m*umpix) ypix(n)+(m*umpix)],'-y')    %% plotting x positions of punctum start+end back onto image
    plot([xt2(n) xt2(n)],[ypix(n)-(m*umpix) ypix(n)+(m*umpix)],'-y')

    %%% showing Ch2 image
    Ch2 = I(j(n)).data; %%% specifying the analysis image as 'Ch2'
    mxfj(n) = max(max(Ch2));
    subplot(2,3,5)
    imshow(Ch2,[0 mxfj(n)])
    title('Ch2 zoom')
    colormap('jet')
    if nargin > 2
        xpix(n) = xpix(n)+(translate(1)*umpix);   %%% offsetting 10um for translation control
        xt1(n) = xt1(n)+(translate(1)*umpix);
        xt2(n) = xt2(n)+(translate(1)*umpix);
        ypix(n) = ypix(n)+(translate(2)*umpix);
    end
    axis([xpix(n)-(m*umpix) xpix(n)+(m*umpix) ypix(n)-(m*umpix) ypix(n)+(m*umpix)]) %%% zoomed-in image of punctum +/-m*1um from centroid
    hold on
    plot([xpix(n)-(m*umpix) xpix(n)+(m*umpix)],[ypix(n)+3 ypix(n)+3],'-y')
    plot([xpix(n)-(m*umpix) xpix(n)+(m*umpix)],[ypix(n)-3 ypix(n)-3],'-y')
    plot([xt1(n) xt1(n)],[ypix(n)-(m*umpix) ypix(n)+(m*umpix)],'-y')    %%% plotting x positions of punctum start+end onto Ch2 image
    plot([xt2(n) xt2(n)],[ypix(n)-(m*umpix) ypix(n)+(m*umpix)],'-y')
    
    %%% Ch2 fluorescence within Ch1-defined punctum
    Ch2_linex{n} = [round(xpix(n)-(m*umpix)):round(xpix(n)+(m*umpix))];
    Ch2_liney{n} = round(ypix(n))*ones(1,length(Ch2_linex{n}));
    if min(Ch2_linex{n})>0 & min(Ch2_liney{n})>0 & max(Ch2_linex{n})<2049 & max(Ch2_liney{n})<2049 %%% check in case translation has gone over the image edge! Returns 'cannot measure' in this case
    for k = 1:length(Ch2_linex{n})
        Ch2_linef{n}(k) = Ch2(round(ypix(n)),Ch2_linex{n}(k)); %%%% matlab uses rows then columns, so this has to be in [y,x] co-ordinates!
        Ch2_linef3{n}(k) = mean([Ch2(round(ypix(n))-1,Ch2_linex{n}(k)) Ch2(round(ypix(n)),Ch2_linex{n}(k)) Ch2(round(ypix(n))+1,Ch2_linex{n}(k))]); %%% mean of 3 pixels centred on punctum y-centre
        Ch2_linef7{n}(k) = mean([Ch2(round(ypix(n))-3,Ch2_linex{n}(k)) Ch2(round(ypix(n))-2,Ch2_linex{n}(k)) Ch2(round(ypix(n))-1,Ch2_linex{n}(k)) Ch2(round(ypix(n)),Ch2_linex{n}(k)) Ch2(round(ypix(n))+1,Ch2_linex{n}(k)) Ch2(round(ypix(n))+2,Ch2_linex{n}(k)) Ch2(round(ypix(n))+3,Ch2_linex{n}(k))]); %%% mean of 7 pixels centred on punctum y-centre
    end
    
    subplot(2,3,6)
    plot(Ch2_linex{n},Ch2_linef7{n},'g-')
    xlabel('pixels')
    ylabel('F')
    hold on

    Ch2_prex{n} = [(xt1(n)-1-(1*umpix)):(xt1(n)-1)];  %%% finding x pixels 2um before Ch1-defined punctum start
    for k = 1:length(Ch2_prex{n})
        Ch2_pref{n}(k) = mean([Ch2(round(ypix(n))-3,Ch2_prex{n}(k)) Ch2(round(ypix(n))-2,Ch2_prex{n}(k)) Ch2(round(ypix(n))-1,Ch2_prex{n}(k)) Ch2(round(ypix(n)),Ch2_prex{n}(k)) Ch2(round(ypix(n))+1,Ch2_prex{n}(k)) Ch2(round(ypix(n))+2,Ch2_prex{n}(k)) Ch2(round(ypix(n))+3,Ch2_prex{n}(k))]); %%% mean of 7 pixels centred on punctum y-centre
    end
        
    Ch2_postx{n} = [(xt2(n)+1):(xt2(n)+(1*umpix))];  %%% finding x pixels 2um before Ch1-defined punctum start
    for k = 1:length(Ch2_postx{n})
        Ch2_postf{n}(k) = mean([Ch2(round(ypix(n))-3,Ch2_postx{n}(k)) Ch2(round(ypix(n))-2,Ch2_postx{n}(k)) Ch2(round(ypix(n))-1,Ch2_postx{n}(k)) Ch2(round(ypix(n)),Ch2_postx{n}(k)) Ch2(round(ypix(n))+1,Ch2_postx{n}(k)) Ch2(round(ypix(n))+2,Ch2_postx{n}(k)) Ch2(round(ypix(n))+3,Ch2_postx{n}(k))]); %%% mean of 7 pixels centred on punctum y-centre
    end
        
    
    Ch2_prefmean(n) = mean(Ch2_pref{n});
    Ch2_postfmean(n) = mean(Ch2_postf{n});
    [Ch2_backg(n) prepost(n)] = min([Ch2_prefmean(n) Ch2_postfmean(n)]);  %% finding background (pre vs post) with lower mean F
    if prepost(n)==1
        plot([Ch2_prex{n}(1) Ch2_prex{n}(end)],[Ch2_backg(n) Ch2_backg(n)],'m-')
        plot(Ch2_prex{n},Ch2_pref{n},'m-')
        Ch2_backg_sd(n) = std([Ch2_pref{n}]); %%% standard deviation of Ch2 background F
        Ch2_threshf(n) = Ch2_backg(n)+(1*Ch2_backg_sd(n));
        plot([Ch2_linex{n}(1) Ch2_linex{n}(end)],[Ch2_threshf(n) Ch2_threshf(n)],'m-')
    elseif prepost(n)==2
        plot([Ch2_postx{n}(1) Ch2_postx{n}(end)],[Ch2_backg(n) Ch2_backg(n)],'m-')
        plot(Ch2_postx{n},Ch2_postf{n},'m-')
        Ch2_backg_sd(n) = std([Ch2_postf{n}]); %%% standard deviation of Ch2 background F
        Ch2_threshf(n) = Ch2_backg(n)+(1*Ch2_backg_sd(n));
        plot([Ch2_linex{n}(1) Ch2_linex{n}(end)],[Ch2_threshf(n) Ch2_threshf(n)],'m-')
    end
    
    punct_i{n} = find((Ch2_linex{n}>=xt1(n))&(Ch2_linex{n}<=xt2(n))); %% finding Ch2 F over all Ch1-defined punctum
    Ch2_punctf{n} = Ch2_linef7{n}(punct_i{n});
    Ch2_punctx{n} = Ch2_linex{n}(punct_i{n});
    [Ch2_punctf_mx(n) Ch2_punctf_mxi(n)] = max(Ch2_punctf{n}); %%% maxF of Ch2 within Ch1-defined punctum
    Ch2_punctf_mxx(n) = Ch2_punctx{n}(Ch2_punctf_mxi(n));
    plot(Ch2_punctf_mxx(n),Ch2_punctf_mx(n),'bo')
    Ch2_maxx{n} = (Ch2_punctf_mxx(n)-(umpix/2)):(Ch2_punctf_mxx(n)+(umpix/2)-1);
    for k = 1:length(Ch2_maxx{n})
        Ch2_maxf{n}(k) = mean([Ch2(round(ypix(n))-3,Ch2_maxx{n}(k)) Ch2(round(ypix(n))-2,Ch2_maxx{n}(k)) Ch2(round(ypix(n))-1,Ch2_maxx{n}(k)) Ch2(round(ypix(n)),Ch2_maxx{n}(k)) Ch2(round(ypix(n))+1,Ch2_maxx{n}(k)) Ch2(round(ypix(n))+2,Ch2_maxx{n}(k)) Ch2(round(ypix(n))+3,Ch2_maxx{n}(k))]); %%% mean of 7 pixels centred on punctum y-centre
    end
    plot(Ch2_maxx{n},Ch2_maxf{n},'b-')
    Ch2_F(n) = mean(Ch2_maxf{n});
    plot([Ch2_maxx{n}(1) Ch2_maxx{n}(end)],[Ch2_F(n) Ch2_F(n)],'b-')
    
    if Ch2_F(n)>Ch2_threshf(n)
        cloc(n) = 1;
        disp('colocalised')
        title(['Ch2 profile: colocalised'])
    else
        cloc(n) = 0;
        disp('not colocalised')
        title(['Ch2 profile: not colocalised'])
    end

    hold off

    pause
    else
        cloc(n) = 999;
        disp('cannot measure')
    end
end

cloc = cloc'

perc_cloc = mean(cloc(find(cloc~=999)))*100





