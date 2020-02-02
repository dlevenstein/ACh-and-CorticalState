analysisfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/AnalysisScripts/AnalysisFigs/BehaviorAChAnalysis';
BehaviorAChAll = GetMatResults(analysisfolder,'BehaviorAChAnalysis');
%genotype = {PupilWhiskAll.genotype};
%[genotypes,~,genotypeidx] = unique(genotype);
%%
BehaviorACh = bz_CollapseStruct(BehaviorAChAll);

%%

pupilphaseonsetGACh = bz_CollapseStruct(BehaviorACh.pupilphaseonsetGACh,3,'mean',true);
pupilphaseGACh = bz_CollapseStruct(BehaviorACh.pupilphaseGACh,3,'mean',true);
pupildpGACh = bz_CollapseStruct(BehaviorACh.pupildpGACh,3,'mean',true);
GAChdist = bz_CollapseStruct(BehaviorACh.GAChdist,3,'mean',true);
PAcoupling = bz_CollapseStruct(BehaviorACh.PAcoupling,4,'mean',true);
%%
HILO = {'lopup','hipup'};
WHNWH = {'Wh','NWh'};
ONOFF = {'WhOn','WhOFF'};
LONGSHORT = {'long','short'};
GAChrange = [0.85 1.2];
%%
%%
figure
for oo = 1:2
    subplot(2,2,oo)
        imagesc(pupilphaseonsetGACh.Xbins,pupilphaseonsetGACh.Ybins,...
            pupilphaseonsetGACh.(ONOFF{oo}).meanZ')
        alpha(single(~isnan(pupilphaseonsetGACh.(ONOFF{oo}).meanZ')))
        hold on
        axis xy
        plot([0 0],ylim(gca),'w--')
        xlabel(['t - aligned to ',(ONOFF{oo})]);ylabel('Pupil Phase')
        ColorbarWithAxis([0.9 1.1],'Mean GACh')
        crameri('berlin','pivot',1)
end
NiceSave('GAChWhisk_aligned',analysisfolder,'')

 
 %%
 figure
for oo = 1:2
    subplot(4,2,oo)
        imagesc( GAChdist.(ONOFF{oo}).all.Xbins,...
            GAChdist.(ONOFF{oo}).all.Ybins,...
            GAChdist.(ONOFF{oo}).all.pYX')
        hold on; axis xy; box off
        colorbar
        plot([0 0],ylim(gca),'k')
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel(['t - ',(ONOFF{oo})]);ylabel('GACh')
   crameri bilbao
   
   for ll = 1:2
        subplot(4,2,oo+ll*2)
            imagesc( GAChdist.(ONOFF{oo}).(LONGSHORT{ll}).Xbins,...
                GAChdist.(ONOFF{oo}).(LONGSHORT{ll}).Ybins,...
                GAChdist.(ONOFF{oo}).(LONGSHORT{ll}).pYX')
            hold on; axis xy; box off
            colorbar
            plot([0 0],ylim(gca),'k')
            %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
            xlabel(['t - ',(ONOFF{oo})]);ylabel('GACh')
       crameri bilbao
   end
end

    subplot(4,2,7)
        imagesc( GAChdist.EMG.Xbins,...
            GAChdist.EMG.Ybins,...
            GAChdist.EMG.pYX')
        hold on; axis xy; box off
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('EMG');ylabel('GACh')
   crameri bilbao
NiceSave('GAChWhisk',analysisfolder,'')
  
  
  %%
cosx = linspace(-pi,pi,100);
cospamp = [0.08 0.8];
figure


subplot(3,2,1)
a = imagesc(pupilphaseGACh.Xbins,pupilphaseGACh.Ybins,pupilphaseGACh.meanZ');
hold on
alpha(a,double(~isnan(pupilphaseGACh.meanZ')))

%imagesc(pupilphaseEMG.Xbins+2*pi,pupilphaseEMG.Ybins,pupilphaseEMG.meanZ')
crameri lapaz
%plot([-pi 3*pi],pupilcycle.detectionparms.pupthresh.*[1 1],'w--')
plot(cosx,(cos(cosx)+1).*cospamp(1)-2,'k')
axis xy
box off
%xlim([-pi 3*pi])
ColorbarWithAxis(GAChrange,'Mean GACh')
%LogScale('c',10)
LogScale('y',10)
xlabel('Pupil Phase');ylabel('Pupil Amplitude')

subplot(3,2,2)
a = imagesc(pupildpGACh.Xbins,pupildpGACh.Ybins,pupildpGACh.meanZ');
hold on
alpha(a,double(~isnan(pupildpGACh.meanZ')))
plot(pupildpGACh.Xbins([1 end]),[0 0],'k--')
crameri lapaz
ColorbarWithAxis(GAChrange,'Mean GACh')
%LogScale('c',10)
axis xy
box off
%LogScale('c',10)
LogScale('x',10)



subplot(3,2,3)
        for pp = 1:2
        imagesc( GAChdist.(HILO{pp}).Xbins+2*pi*(pp-1),...
            GAChdist.(HILO{pp}).Ybins,...
            GAChdist.(HILO{pp}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-2,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        colorbar
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel('GACh')
        crameri bilbao
        
subplot(3,2,4)

        imagesc( GAChdist.pup.Xbins,...
            GAChdist.pup.Ybins,...
            GAChdist.pup.pYX')
        hold on; axis xy; box off
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('GACh')
   crameri bilbao
NiceSave('GAChPupil',analysisfolder,'')

   %% Figure: whisk/nonwhisk
   figure
for ww = 1:2
 subplot(3,2,ww)

        imagesc( GAChdist.pup.(WHNWH{ww}).Xbins,...
            GAChdist.pup.(WHNWH{ww}).Ybins,...
            GAChdist.pup.(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        colorbar
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        xlabel('Pupil Size');ylabel('GACh')
        title((WHNWH{ww}))
   crameri bilbao
   
subplot(3,2,2+ww)
        for pp = 1:2
        imagesc( GAChdist.(HILO{pp}).(WHNWH{ww}).Xbins+2*pi*(pp-1),...
            GAChdist.(HILO{pp}).(WHNWH{ww}).Ybins,...
            GAChdist.(HILO{pp}).(WHNWH{ww}).pYX')
        hold on; axis xy; box off
        plot(cosx+2*pi*(pp-1),(cos(cosx)+1).*cospamp(pp)-2,'k')
        end   
        %ColorbarWithAxis([-2.4 -1.2],'Mean PSS')
        colorbar
        xlim([-pi 3*pi])
        xlabel('Pupil Phase');ylabel('GACh')
        crameri bilbao
        
end
   
NiceSave('GAChPupil_WhNWh',analysisfolder,'')
%%
%%
figure
for oo = PAcoupling.orders
    subplot(2,3,oo)
imagesc(log10(PAcoupling.lowerbounds),log10(PAcoupling.upperbounds),PAcoupling.coupling(:,:,oo)')
hold on
axis xy
%plot(log10(pupilcycle.detectionparms.filterparms.passband(1)),log10(pupilcycle.detectionparms.filterparms.passband(2)),'r+')
LogScale('xy',10)
clim([0.01 0.023])
xlabel('Lower Bound (Hz)');ylabel('Upper Bound (Hz)')
colorbar
title(num2str(oo))


    subplot(2,3,oo+3)
imagesc(log10(PAcoupling.lowerbounds),log10(PAcoupling.upperbounds),PAcoupling.corr(:,:,oo)')
hold on
axis xy
%plot(log10(pupilcycle.detectionparms.filterparms.passband(1)),log10(pupilcycle.detectionparms.filterparms.passband(2)),'r+')
LogScale('xy',10)
%clim([0.1 0.4])
xlabel('Lower Bound (Hz)');ylabel('Upper Bound (Hz)')
colorbar
caxis([0.1 0.6])
title(num2str(oo))
end

NiceSave('FilterCompare',analysisfolder,'')