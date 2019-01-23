function [  ] = LaminarCSDAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = pwd;
figfolder = '/mnt/data1/Dropbox/research/Current Projects/S1State/AnalysisScripts/figures/UPDOWNandPupilAnalysis';

%%
baseName = bz_BasenameFromBasepath(basePath);
recparms = bz_getSessionInfo(basePath,'noPrompts',true);

%% Detect Slow Waves
%CTXChans = recparms.SpkGrps.Channels(23:46);
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');

%%
pupildilation = bz_LoadBehavior(basePath,'pupildiameter');



end

