function xdot = SeasonLV_TypeII(tt, xt, SpParam, IntParam, IsStep, IsReproSeasonal)
	% SeasonLV_TypeII() returns the growth rates of species involved in a seasonal food web described by IntParam
	% the functional response is of type II
	% the parameters relative to the intra-specific processes are described in SpParam
	% the mortality of predators is density dependent (and tuned by parameter rho)
	% the type of forcing signal is described with logical IsStep
	% prey reproduction can be set as seasonal with logical IsReproSeasonal

	% IntParam is a table listing the interactions
	if (isempty(IntParam.LowerTaxon)) | (isempty(IntParam.UpperTaxon)) | (isempty(IntParam.E)) | (isempty(IntParam.G)) | (isempty(IntParam.H)) | (isempty(IntParam.SeasonDom))
		error('IntParam must have at least 6 columns, namely `LowerTaxon`, `UpperTaxon`, `E`, `G`, `H`, and `SeasonDom`.');
	end
	% SpParam is a table listing each species parameter
	if (isempty(SpParam.Taxon)) | (isempty(SpParam.TrophicLevel)) | (isempty(SpParam.R)) | (isempty(SpParam.M)) | (isempty(SpParam.BETA) | (isempty(SpParam.IsCst)) | (isempty(SpParam.DD)))
		error('SpParam must have at least 7 columns, namely `Taxon`, `TrophicLevel`, `M`, `R`, `BETA`, `IsCst`, and `DD`.');
	end
	% IsStep must be a logical variable
	if ~islogical(IsStep)
		error('IsStep must be a logical variable.')
	end
	% IsReproSeasonal must be a logical variable
	if ~islogical(IsReproSeasonal)
		error('IsReproSeasonal must be a logical variable.')
	end

	global EXTINCT_THRS CONV_EFF


	if IsStep
		ForcingSignalS = 1;
		ForcingSignalW = -1;

		Year = floor(tt); % year
		ttprime = tt-Year; % time of the year

		if ttprime >= 0.5
			% if winter
			ForcingSignalS = -1;
			ForcingSignalW = 1;
		end
	else
		ForcingSignalS = sin(2*pi*tt);
		ForcingSignalW = sin(2*pi*tt + pi);
	end


	STot = length(xt); % total number of species
	SPrey = sum(ismember(SpParam.TrophicLevel, 'Prey'));
	SPred = sum(ismember(SpParam.TrophicLevel, 'Predator'));
	xdot = zeros(STot, 1);

	xt(xt < EXTINCT_THRS) = 0;
	for sp = 1:STot
		if (xt(sp) >= EXTINCT_THRS) & (~SpParam.IsCst(sp))
			if string(SpParam.TrophicLevel(sp)) == 'Prey'

				% calculate growth rate caused by reproduction and intra-specific competition
				if IsReproSeasonal
%					xdot(sp) = xdot(sp) + (SpParam.R(sp) - SpParam.BETA(sp)*xt(sp)) * (1 + ForcingSignalS); % assuming reproduction happens during summer
					xdot(sp) = xdot(sp) + SpParam.R(sp) * (1 + ForcingSignalS) - SpParam.BETA(sp)*xt(sp); % assuming reproduction happens during summer
				else
					xdot(sp) = xdot(sp) + SpParam.R(sp) - SpParam.BETA(sp)*xt(sp);
				end

				int_sp = find(ismember(IntParam.LowerTaxon, SpParam.Taxon(sp))); % list of species in which sp is involved as a prey
				for int = int_sp.'
					pred = IntParam.UpperTaxon(int); % predator's name
					predIndexSp = find(ismember(SpParam.Taxon, pred)); % index of pred in SpParam


				% calculate sum of discovered prey biomass sum_j a_ij(t) * x_j
					listOfPrey = IntParam.LowerTaxon(ismember(IntParam.UpperTaxon, pred));
					int_pred = find(ismember(IntParam.UpperTaxon, pred)); % index of pred's interactions
					sumDiscoveredPreysB = 0;
					for int = int_sp.'
						prey = IntParam.LowerTaxon(int); % preys' names
						preyIndexSp = find(ismember(SpParam.Taxon, prey)); % index of preys in SpParam

						preyBiomass = 1; % prey biomass (default is 1)
						if (~SpParam.IsCst(preyIndexSp))
							preyBiomass = xt(preyIndexSp);
						end

						if string(IntParam.SeasonDom(int)) == 'S'
							sumDiscoveredPreysB = sumDiscoveredPreysB + IntParam.G(int) * (1 + IntParam.E(int) * ForcingSignalS) * preyBiomass;
						else
							sumDiscoveredPreysB = sumDiscoveredPreysB + IntParam.G(int) * (1 + IntParam.E(int) * ForcingSignalW) * preyBiomass;
						end
						if isnan(sumDiscoveredPreysB)
							disp([sp prey xt(preyIndexSp) int]);
							error("Problem in the sum of discovered prey biomass.");
						end
					end
					
					% calculate a_ik / (1 + h_i * sum_j(a_ij(t) * x_j(t)))
					FR = IntParam.G(int) / (1 + IntParam.H(int) * sumDiscoveredPreysB);

					% FR * xt(predIndexSp) * (1 + IntParam.E(int)*ForcingSignalS) / SpParam.BodyMass(predIndexSp) = a_ki(t) * x_i(t) / (1 + h_i * sum_j(a_ij(t) * xj(t)))
					if string(IntParam.SeasonDom(int)) == 'S'
						xdot(sp) = xdot(sp) - FR * xt(predIndexSp) * (1 + IntParam.E(int)*ForcingSignalS) / SpParam.BodyMass(predIndexSp); % if summer interaction
					else
						xdot(sp) = xdot(sp) - FR * xt(predIndexSp) * (1 + IntParam.E(int)*ForcingSignalW) / SpParam.BodyMass(predIndexSp); % if winter interaction
					end
				end
			else
				xdot(sp) = xdot(sp) - SpParam.M(sp) - SpParam.DD(sp)*xt(sp);
				int_sp = find(ismember(IntParam.UpperTaxon, SpParam.Taxon(sp)));

				% calculate sum of discovered prey biomass sum_j a_ij(t) * x_j
				listOfPrey = IntParam.LowerTaxon(int_sp);
				indexPrey = find(ismember(SpParam.Taxon, listOfPrey));
				sumDiscoveredPreysB = 0;
				for int = int_sp.'
					prey = IntParam.LowerTaxon(int); % preys' names
					preyIndexSp = find(ismember(SpParam.Taxon, prey)); % index of preys in SpParam

					preyBiomass = 1; % prey biomass (default is 1)
					if (~SpParam.IsCst(preyIndexSp))
						preyBiomass = xt(preyIndexSp);
					end

					if string(IntParam.SeasonDom(int)) == 'S'
						sumDiscoveredPreysB = sumDiscoveredPreysB + IntParam.G(int) * (1 + IntParam.E(int) * ForcingSignalS) * preyBiomass;
					else
						sumDiscoveredPreysB = sumDiscoveredPreysB + IntParam.G(int) * (1 + IntParam.E(int) * ForcingSignalW) * preyBiomass;
					end
					if isnan(sumDiscoveredPreysB)
						disp([sp prey xt(preyIndexSp) int]);
						error("Problem in the sum of discovered prey biomass.");
					end
				end

				for int = int_sp.'
					prey = IntParam.LowerTaxon(int); % preys' names
					preyIndexSp = find(ismember(SpParam.Taxon, prey)); % index of preys in SpParam

					% calculate a / (1 + h * SumPreysB * a) * 1 / Mpred
					FR = IntParam.G(int) / (1 + IntParam.H(int) * sumDiscoveredPreysB);

					if string(IntParam.SeasonDom(int)) == 'S'
						if (~SpParam.IsCst(preyIndexSp))
							xdot(sp) = xdot(sp) + CONV_EFF * FR * xt(preyIndexSp) * (1 + IntParam.E(int)*ForcingSignalS) / SpParam.BodyMass(sp);
						else
							xdot(sp) = xdot(sp) + CONV_EFF * FR * (1 + IntParam.E(int)*ForcingSignalS) / SpParam.BodyMass(sp);
						end
						if isnan(xdot(sp))
							disp([xdot(sp) prey xt(preyIndexSp) int IntParam.E(int) IntParam.G(int) IntParam.H(int) sumDiscoveredPreysB FR sp]);
							error("truc");
						end
					else
						if (~SpParam.IsCst(preyIndexSp))
							xdot(sp) = xdot(sp) + CONV_EFF * FR * xt(preyIndexSp) * (1 + IntParam.E(int)*ForcingSignalW) / SpParam.BodyMass(sp);
						else
							xdot(sp) = xdot(sp) + CONV_EFF * FR * (1 + IntParam.E(int)*ForcingSignalW) / SpParam.BodyMass(sp);
						end
						if isnan(xdot(sp))
							disp([xdot(sp) prey xt(preyIndexSp) IntParam.E(int) IntParam.G(int) IntParam.H(int) sumDiscoveredPreysB FR sp]);
							error("truc");
						end
					end
				end
			end
			xdot(sp) = xdot(sp)*xt(sp);
		end
	end
end
