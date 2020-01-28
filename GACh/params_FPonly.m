function [ params ] = params_FPonly
%% Process Parameters for Photometry Analysis
%Created By: Pratik Mistry
%Created On: 31 January 2019
%Edited On: 12 June 2019
%
%Description: This is a script with different variables whose values can be
%adjusted depending on the photometry signal that is being processed. The
%name of this file can and should be changed depending on the method and
%GECI used.
%
%
%% General Parameters
params.dsRate = 0; % Downsampling rate if you want to downsample the signal
%This dsRate will also be applied to all signals during the analysis
%pipeline

%% Demodulation Parameters
%Adjust the demodStatus variable to "1" if you need to demodulate a signal
%from a lock-in amplifier or "0" if it's a normal photometry recording
params.FP.demodStatus = 0; 

%% Filter Parameters
%params.filtType = 'lowpass'; % Filter type: 'lowpass' or 'highpass' --
%Temporarily removed. We more than likely won't need to highpass filter our
%signals
params.FP.lpCut = 10; % Cut-off frequency for filter
params.FP.filtOrder = 10; % Order of the filter

%% Baseline Parameters
params.FP.interpType = 'linear'; % 'linear' 'spline' 
params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'
params.FP.winSize = 20; % Window size for baselining in seconds
params.FP.winOv = 1; %Window overlap size in seconds
params.FP.basePrc = 10; % Percentile value from 1 - 100 to use when finding baseline points