
figfolder = '/Users/dlevenstein/Project Repos/ACh-and-CorticalState/Modeling/ModelingFigs';


%% Load I Bif'n

addpath('/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/Modeling/simresults');

BifI = 5;

Ibifn = importdata(['Ibifn_w',num2str(BifI),'.dat']);
LCbifn = importdata(['Ibifn_w',num2str(BifI),'_LC.dat']);
    
stableline = [Ibifn(:,1) Ibifn(:,3)];
stableline(Ibifn(:,4)==2,:) =nan;

unstableline = [Ibifn(:,1) Ibifn(:,3)];
unstableline(Ibifn(:,4)==1,:) =nan;

numpoints = 50;
LCline = [LCbifn(22:2:numpoints,1),LCbifn(22:2:numpoints,2),LCbifn(22:2:numpoints,3)];

%%
%Double noise UP/DOWN model
simtime = 50000;
dt = 1;

parms.N_neurons = 1;
parms.I_in = -1.65;
parms.W = 5;
parms.beta = 1;
parms.tau_r = 1;
parms.tau_a = 25;
parms.A0 = 0.5;
parms.Ak = 15;
parms.noiseamp = [0.1 0.1];
parms.noisefreq = [1 0.005];
[ T, Y_sol,Inoise,Ipulse ] = WCadapt_run(simtime,dt,parms);


%Calculate DWELL times
[thresh,cross,~,diptest] = BimodalThresh(Y_sol(:,1),'Schmidt');

if isempty(cross.upints) || length(cross.upints) <=2
    dwell.UP = nan; dwell.DOWN = nan;
else
    dwell.UP = cross.upints(:,2)-cross.upints(:,1);
    dwell.DOWN = cross.downints(:,2)-cross.downints(:,1);
end

%%
UPDOWN = {'UP','DOWN'};
UDcolor = {'r','b'};

%%
winsize = 5000;
xwin = 100+[0 winsize];
figure
    subplot(4,3,1:2)
        hold on
        plot(stableline(:,1),stableline(:,2),'k','Linewidth',2)
        plot(unstableline(:,1),unstableline(:,2),'k:','Linewidth',1)
        plot(LCline(:,1),LCline(:,2),'ko','markersize',3)
        plot(LCline(:,1),LCline(:,3),'ko','markersize',3)
        plot(parms.I_in,1,'k+')
        plot(parms.I_in+parms.noiseamp(2)*[-1 1],[1 1],'k')
        xlim([-2.5 -1])
        xlabel('Drive');ylabel('R_s_s')
    
    subplot(4,1,2)
        plot(T,Y_sol(:,1),'k','linewidth',1)
        box off
        ylabel('Rate')
        xlim(xwin)
    
    subplot(4,1,3)
        plot(T,sum(Inoise,3),'color',[0.5 0.5 0.5],'linewidth',0.5)
        hold on
        plot(T,Inoise(:,:,2),'k','linewidth',1)
        box off
        ylabel('Noise')
        xlim(xwin)
    
    for uu = 1:2
    subplot(4,4,12+uu)
        plot(log10(dwell.(UPDOWN{uu})(1:end-1)),log10(dwell.(UPDOWN{uu})(2:end)),'.','color',UDcolor{uu})
        
        xlim([1 3]);ylim([1 3])
        LogScale('xy',10)
    end
    
    bins = linspace(1,3,30);
    subplot(4,3,3)
    hold on
        for uu = 2:-1:1
            counts = hist(log10(dwell.(UPDOWN{uu})),bins);
            bar(bins,counts,'FaceColor',UDcolor{uu})
        end
        LogScale('x',10)
        
        NiceSave('SlowFluct',figfolder,'WCadapt')
