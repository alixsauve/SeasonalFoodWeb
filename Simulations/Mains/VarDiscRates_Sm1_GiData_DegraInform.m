% Scenario "Degrading information":
% We alter the densities of the best-studied taxonomic groups (namely mammals and birds by increasing or decreasing their initial densities by 10% or 50%.

clc
clear all
close all

DDName = 'DDData';

rng(2019); % set random seed

global EXTINCT_THRS CONV_EFF
CONV_EFF = 0.1;
EXTINCT_THRS = 1e-6;

% defining length of simulations, and time steps
TFinal = 100; % 100 years are simulated
TStep = 1/52; % we want the densities every weeks
TSpan = 0:TStep:TFinal; % vector of time steps we want as outputs of simulations

% add repertories to the working directory
addpath('../Functions');

% load reference parameters
% parameter values for species intra-specific processes
ReadOptions = detectImportOptions('../../Parameterisation/OutputTables/YearPopParam.csv');
ReadOptions = setvartype(ReadOptions, {'M', 'DDData', 'DDModel'}, 'double');
PopParamRef = readtable('../../Parameterisation/OutputTables/YearPopParam.csv', ReadOptions);
PopParamRef.Properties.VariableNames{DDName} = 'DD'; % rename the DD column to be used with 'DD' (i.e., the variable name in growth rate functions).

% parameters of predatory links, type I functional response
ReadOptions = detectImportOptions('../../Parameterisation/OutputTables/YearIntParam_TypeI.csv');
IntParamRefI = readtable('../../Parameterisation/OutputTables/YearIntParam_TypeI.csv', ReadOptions);

% initial biomass densities (mean values observed on the field)
ReadOptions = detectImportOptions('../../Parameterisation/OutputTables/SpDensBiomass.csv');
ICRef = readtable('../../Parameterisation/OutputTables/SpDensBiomass.csv', ReadOptions);
ICRef.InitBiomass_gha = ICRef.BodyMass_g.*ICRef.InitDensity_Nha;

% information on species
NSp = size(PopParamRef, 1); % number of species
PreysIndex = ismember(PopParamRef.TrophicLevel, 'Prey');
PredsIndex = ismember(PopParamRef.TrophicLevel, 'Predator');

% number of replicates per case
NRep = 3;

% options for the numerical integration
Options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'NonNegative', [1:NSp]);

GRFxn = @SeasonLV_TypeI; % growth rates function
IntParam = IntParamRefI; % predation parameters 

% indices of taxa for which we change initial densities
% indices of ungulates
UnguIndices = ismember(ICRef.Clade, 'Ungulate');
UnguList = ICRef.Taxon(UnguIndices);
% indices of rodents
RodeIndices = ismember(ICRef.Clade, 'Rodent');
RodeList = ICRef.Taxon(RodeIndices);
% indices of birds
BirdIndices = logical(ismember(ICRef.CladeSimpler, 'Bird') .* PreysIndex);
BirdList = ICRef.Taxon(BirdIndices);
% indices of species from these three groups
AggIndices = logical(BirdIndices + RodeIndices + UnguIndices);
AggList = ICRef.Taxon(AggIndices);

SimuTime = []; % an empty vector to display simulation times

% prey reproduce only during one season
for Group = {'Agg' 'Bird' 'Rode' 'Ungu'}
	GroupIndices = eval(string(strcat(Group, 'Indices'))); % get species indices in species list
	for Rep = 1:NRep
		Sign = (-1) .^ randi([1 2], size(ICRef, 1), 1); % pick whether to increase or decrease initial densities
		for Fact = [0.1 0.5]
			disp(strcat('Degrading information; GiData; Replicate #', num2str(Rep), '; FR type ', func2str(GRFxn), '; Deviation: ', num2str(Fact), '; Targetted group: ', string(Group)));

			GroupIndices = eval(string(strcat(Group, 'Indices'))); % get species indices in species list
			GroupSpList = eval(string(strcat(Group, 'List'))); % get the species list (names)
			GroupSize = sum(GroupIndices); % number of species in that group

			% change discovery rate on each species of this group in respect with the way densities are altered
			for Sp = 1:GroupSize
				SpIntIndices = ismember(IntParam.LowerTaxon, GroupSpList{Sp}); % indices of species interactions
				IntParam.G(SpIntIndices) = IntParam.G(SpIntIndices) / (1 + Sign(Sp) * Fact); % modify the corresponding discovery rates
			end

			% now running the simulations
			tic;
			[T, Y] = ode45(@(t, y) GRFxn(t, y, PopParamRef, IntParam, false, true), TSpan, ICRef.InitBiomass_gha, Options);
			SimuTime = [SimuTime toc]


			% calculate the persistence
			sum(Y(end, :) >= EXTINCT_THRS)/NSp
		
			% save the time series
			dlmwrite(char(strcat('../Outputs/TS_VarDiscRates_Sm1_', DDName, '_DegraInform_Rep', num2str(Rep), '_Dev', num2str(Fact), '_', string(Group), '.csv')), [T Y], 'delimiter', ',');
			clear T Y
		end
	end
end
