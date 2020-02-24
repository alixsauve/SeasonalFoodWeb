% This script runs the dynamics of the seasonal food web of the Bialowieza forest following two types of predator functional responses.

% All are ran with a sinusoidal forcing signal.
% For all of these simulations, the predator density-dependent mortality rates are based on the maximum prey intake estimated for each predator (DD = DDModel).

clc
clear all
close all

DDName = 'DDModel';

rng(2019); % set random seed

global EXTINCT_THRS CONV_EFF
CONV_EFF = 0.1;
EXTINCT_THRS = 1e-6;

% defining length of simulations, and time steps
TFinal = 100; % 100 years are simulated
TStep = 1/52; % we want the densities every weeks
TSpan = 0:TStep:TFinal; % vector of time steps we want as outputs of simulations

% add repertories to the working directory
addpath('../Functions'); % path to functions

% load reference parameters
% parameter values for species intra-specific processes
ReadOptions = detectImportOptions('../../Parameterisation/OutputTables/YearPopParam.csv');
ReadOptions = setvartype(ReadOptions, {'M', 'DDModel', 'DDData'}, 'double');
PopParamRef = readtable('../../Parameterisation/OutputTables/YearPopParam.csv', ReadOptions);
PopParamRef.Properties.VariableNames{DDName} = 'DD'; % rename the DD column to be used with 'DD' (i.e., the variable name in growth rate functions).

% parameters of predatory links, type I functional response
ReadOptions = detectImportOptions('../../Parameterisation/OutputTables/YearIntParam_TypeI.csv');
IntParamRefI = readtable('../../Parameterisation/OutputTables/YearIntParam_TypeI.csv', ReadOptions);

% parameters of predatory links, type II functional response
ReadOptions = detectImportOptions('../../Parameterisation/OutputTables/YearIntParam_TypeII.csv');
IntParamRefII = readtable('../../Parameterisation/OutputTables/YearIntParam_TypeII.csv', ReadOptions);

% initial biomass densities (mean values observed on the field)
ReadOptions = detectImportOptions('../../Parameterisation/OutputTables/SpDensBiomass.csv');
ICRef = readtable('../../Parameterisation/OutputTables/SpDensBiomass.csv', ReadOptions);
ICRef.InitBiomass_gha = ICRef.BodyMass_g.*ICRef.InitDensity_Nha;

% information on species
NSp = size(PopParamRef, 1); % number of species
PreysIndex = ismember(PopParamRef.TrophicLevel, 'Prey');
PredsIndex = ismember(PopParamRef.TrophicLevel, 'Predator');

% options for the numerical integration
Options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'NonNegative', [1:NSp]);
SimuTime = []; % vector of simulation time

% NB: prey reproduce only during one season

% functional response type I
GRFxn = @SeasonLV_TypeI;
disp('Simulating with a type I functional response and a smooth signal, prey reproduction is seasonal...');
[TI, YI] = ode45(@(t, y) GRFxn(t, y, PopParamRef, IntParamRefI, false, true), TSpan, ICRef.InitBiomass_gha, Options);
sum(YI(end, :) >= EXTINCT_THRS)/NSp
dlmwrite(char(strcat('../Outputs/TS_Sm1_', DDName, '.csv')), [TI YI], 'delimiter', ',');
clear TI YI

% functional response type II
GRFxn = @SeasonLV_TypeII;
disp('Simulating with a type II functional response and a smooth signal, prey reproduction is seasonal...');
[TII, YII] = ode45(@(t, y) GRFxn(t, y, PopParamRef, IntParamRefII, false, true), TSpan, ICRef.InitBiomass_gha, Options);
sum(YII(end, :) >= EXTINCT_THRS)/NSp
dlmwrite(char(strcat('../Outputs/TS_Sm2_', DDName, '.csv')), [TII YII], 'delimiter', ',');
clear TII YII
