%%
spathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180708_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180710_WT_EM1M3\Spont1'];

%%
for i = 1:size(spathdatabase,1)
    
    %     basePath = pwd;
    basePath = spathdatabase(i,:);
    cd(basePath);
    
    [baseFolder,baseName] = fileparts(basePath);
    savefile = fullfile(basePath,[baseName,'.LaminarPowerSpec.lfp.mat']);
    figfolder = fullfile(basePath,'AnalysisFigures');
    
    savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
    load(savfile,'sessionInfo');
    badchannels = sessionInfo.badchannels;
    usechannels = sessionInfo.AnatGrps.Channels;
    usechannels(ismember(usechannels,badchannels))=[];
    
    lfp = bz_GetLFP('all','noPrompts',true);
    
    [~,chanidx] = ismember(usechannels,lfp.channels);
    
    clear lfp
    
    %%
    lamwavespec = [];
    lammua = [];
    for ii = 1:length(usechannels)
        [f, pspec, MUA] = bz_MUAGammafromDat(basePath,'channels',chanidx(ii));
        lamwavespec = cat(2,lamwavespec,pspec);
        lammua = cat(2,lammua,[nanmean(MUA); nanmedian(MUA); max(MUA); min(MUA)]);
    end
    
    %     figure;
    %     imagesc(log2(f),[1:size(lamwavespec,2)],lamwavespec')
    %     axis xy
    %     LogScale('x',2)
    %     colormap('jet'); caxis([0 50]);
    %     c = colorbar;
    %     c.Label.String = 'power (dB)';
    %     set(gca,'YDir','reverse');
    %     xlabel('frequency (Hz)'); ylabel('channel');
    
    %% Saving laminar pspectrum
    lampspectrum.data = lamwavespec;
    lampspectrum.muapow = lammua;
    lampspectrum.freq = f;
    
    save(savefile,'lampspectrum')
    
    %% Smoothening
%     smlamwavespec = zeros(size(lamwavespec,1),size(lamwavespec,2));
%     temp_sm = 100; spat_sm = 5;
%     for ch = 1:size(lamwavespec,2)
%         smlamwavespec(:,ch) = smooth(lamwavespec(:,ch),temp_sm,'sgolay');
%     end
%     
%     for t = 1:size(lamwavespec,1)
%         smlamwavespec(t,:) = smooth(smlamwavespec(t,:),spat_sm,'lowess');
%     end

    %% z-scoring across freq
%     zsmlamwavespec = zeros(size(lamwavespec,1),size(lamwavespec,2));
%     for n = 1:size(lamwavespec,1)
%         zsmlamwavespec(n,:) = (smlamwavespec(n,:)-mean(smlamwavespec(n,:)))./std(smlamwavespec(n,:));
%     end
    
    %% Saving figure
%     figure;
%     subplot(2,6,1); hold on;
%     stairs((lammua(1,:)-lammua(1,1))./max(lammua(1,:)-lammua(1,1)),[1:size(lamwavespec,2)],'Color','k');
%     set(gca,'YDir','reverse');
%     xlim([0 1]); ylim([1 size(lamwavespec,2)]);
%     xlabel('norm. MUA power (0.5-5 KHz)'); ylabel('channel');
%     set(gca,'Ytick',[10 20 30 40 50 60]);
%     set(gca,'Yticklabel',{'10','20','30','40','50','60'});
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     box on;
%     
%     subplot(2,6,2); hold on;
%     tempidx = find(f > 100 & f < 200);
%     muapow1 = mean(smlamwavespec(tempidx,:),1);
%     tempidx = find(f > 300 & f < 500);
%     muapow2 = mean(smlamwavespec(tempidx,:),1);
%     tempidx = find(f > 600 & f < 900);
%     muapow3 = mean(smlamwavespec(tempidx,:),1);
%     tempidx = find(f > 1000 & f < 5000);
%     muapow4 = mean(smlamwavespec(tempidx,:),1);
%     stairs(muapow1./max(muapow1),[1:size(lamwavespec,2)],'Color','r');
%     stairs(muapow2./max(muapow2),[1:size(lamwavespec,2)],'Color','g');
%     stairs(muapow3./max(muapow3),[1:size(lamwavespec,2)],'Color','b');
%     stairs(muapow4./max(muapow4),[1:size(lamwavespec,2)],'Color','k');
%     
%     set(gca,'YDir','reverse');
%     xlim([0 1]); ylim([1 size(lamwavespec,2)]);
%     legend({'100-200 Hz','300-500 Hz','600-900 Hz','1-5 KHz'},'Location','southwest');
%     xlabel('norm. power'); 
%     set(gca,'Ytick',[10 20 30 40 50 60]);
%     set(gca,'Yticklabel',{});
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     box on;
%     
%     subplot(2,6,3:4);
%     imagesc(log2(f),[1:size(lamwavespec,2)],smlamwavespec')
%     axis xy
%     LogScale('x',2)
%     colormap('jet'); caxis([0 max(max(smlamwavespec))-30]);
%     c = colorbar;
%     c.Label.String = 'power (dB)';
%     set(gca,'YDir','reverse');
%     set(gca,'Xtick',[log2([1 5 10 25])]);
%     set(gca,'Xticklabel',{'1','5','10','25'});
%     xlim([log2([1 25])]);
%     xtickangle(45);
%     xlabel('frequency (Hz)');
%     set(gca,'Ytick',[10 20 30 40 50 60]);
%     set(gca,'Yticklabel',{});
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     ylim([1 size(lamwavespec,2)]);
%     
%     subplot(2,6,5:6);
%     imagesc(log2(f),[1:size(lamwavespec,2)],smlamwavespec')
%     axis xy
%     LogScale('x',2)
%     colormap('jet'); caxis([0 max(max(smlamwavespec))-30]);
%     c = colorbar;
%     c.Label.String = 'power (dB)';
%     set(gca,'YDir','reverse');
%     set(gca,'Xtick',[log2([1 50 100 250 500 1000 5000])]);
%     set(gca,'Xticklabel',{'1','50','100','250','500','1000','5000'});
%     xtickangle(45);
%     xlabel('frequency (Hz)');
%     set(gca,'Ytick',[10 20 30 40 50 60]);
%     set(gca,'Yticklabel',{});
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     ylim([1 size(lamwavespec,2)]);
%     
%     subplot(2,6,8);hold on;
%     tempidx = find(f > 100 & f < 200);
%     muapow1 = mean(zsmlamwavespec(tempidx,:),1);
%     tempidx = find(f > 300 & f < 500);
%     muapow2 = mean(zsmlamwavespec(tempidx,:),1);
%     tempidx = find(f > 600 & f < 900);
%     muapow3 = mean(zsmlamwavespec(tempidx,:),1);
%     tempidx = find(f > 1000 & f < 5000);
%     muapow4 = mean(zsmlamwavespec(tempidx,:),1);
%     stairs(muapow1./max(muapow1),[1:size(lamwavespec,2)],'Color','r');
%     stairs(muapow2./max(muapow2),[1:size(lamwavespec,2)],'Color','g');
%     stairs(muapow3./max(muapow3),[1:size(lamwavespec,2)],'Color','b');
%     stairs(muapow4./max(muapow4),[1:size(lamwavespec,2)],'Color','k');
%     
%     set(gca,'YDir','reverse');
%     xlim([0 1]); ylim([1 size(lamwavespec,2)]);
%     legend({'100-200 Hz','300-500 Hz','600-900 Hz','1-5 KHz'},'Location','southwest');
%     xlabel('norm. power');
%     set(gca,'Ytick',[10 20 30 40 50 60]);
%     set(gca,'Yticklabel',{'10','20','30','40','50','60'});
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     box on;
%     
%     subplot(2,6,9:10);
%     imagesc(log2(f),[1:size(lamwavespec,2)],zsmlamwavespec')
%     axis xy
%     LogScale('x',2)
%     colormap('jet'); caxis([(max(max(zsmlamwavespec)))*-1 max(max(zsmlamwavespec))]);
%     c = colorbar;
%     c.Label.String = 'power (z-score)';
%     set(gca,'YDir','reverse');
%     set(gca,'Xtick',[log2([1 5 10 25])]);
%     set(gca,'Xticklabel',{'1','5','10','25'});
%     xlim([log2([1 25])]);
%     xtickangle(45);
%     xlabel('frequency (Hz)');
%     set(gca,'Ytick',[10 20 30 40 50 60]);
%     set(gca,'Yticklabel',{});
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     ylim([1 size(lamwavespec,2)]);
%     
%     subplot(2,6,11:12);
%     imagesc(log2(f),[1:size(lamwavespec,2)],zsmlamwavespec')
%     axis xy
%     LogScale('x',2)
%     colormap('jet'); caxis([(max(max(zsmlamwavespec)))*-1 max(max(zsmlamwavespec))]);
%     c = colorbar;
%     c.Label.String = 'power (z-score)';
%     set(gca,'YDir','reverse');
%     set(gca,'Xtick',[log2([1 50 100 250 500 1000 5000])]);
%     set(gca,'Xticklabel',{'1','50','100','250','500','1000','5000'});
%     xtickangle(45);
%     xlabel('frequency (Hz)');
%     set(gca,'Ytick',[10 20 30 40 50 60]);
%     set(gca,'Yticklabel',{});
%     set(gca,'YGrid','on','Layer','top','GridColor',[0 0 0]);
%     ylim([1 size(lamwavespec,2)]);
    
    %%
%     NiceSave('LaminarPowerSpec',figfolder,baseName)
    
%     close all;
    
end