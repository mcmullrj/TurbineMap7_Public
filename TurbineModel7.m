classdef TurbineModel7
    %build model from turbine mapping data
    %estimate turbine map for scaled geometry
    %this version only works for fixed geometry
    %created by Bob McMullen
    
    properties
        TurbineMapData
        TurbineModel
    end
    
    properties(Constant)
        ConvergenceEfficiencyError = 0.0001;
    end
    
    methods
        function obj = TurbineModel7
            obj.TurbineMapData = [];
            obj.TurbineModel = [];
        end
        
        function obj = TurbineMapDataPrep(obj,PathTurbineData,DataType)
            %% read turbine map data
            %older Matlab versions
            %TurbineGeometryInputTable = readtable(PathTurbineData,'Sheet','Flow Area Geometry');
            %newer Matlab versions
            TurbineGeometryInputTable = readtable(PathTurbineData,'Sheet','Flow Area Geometry','NumHeaderLines',1);
            TurbineMapPrep = readmatrix(PathTurbineData,'Sheet','Turbine Map Data');
            TurbineGeometryInput = TurbineGeometryInputTable{:,3};
            TurbineDataSort = sortrows(TurbineMapPrep,1); %sort data by vane position
            %
            %
            %% geometry
            RadiusRotor = TurbineGeometryInput(2)/1000/2; %m
            TrimWheel = TurbineGeometryInput(4);
            RadiusTipRotorOutlet = (RadiusRotor^2*TrimWheel/100)^0.5; %m
            %design variables initial conditions
            AreaRadiusRatioVoluteThroat = TurbineGeometryInput(1)/1000; %m
            WidthRotorInlet = TurbineGeometryInput(3)/1000; %m
            if TurbineGeometryInput(5) > 0
                RadiusNoseRotorOutlet = TurbineGeometryInput(5)/2/1000; %m
            else
                RadiusNoseRotorOutlet = 0.25*RadiusRotor; %m
            end
            AngleBladeRotorOutlet = 55/360*2*pi; %rad
            GeometryTurbineData = [RadiusRotor,RadiusTipRotorOutlet,AreaRadiusRatioVoluteThroat,RadiusNoseRotorOutlet,AngleBladeRotorOutlet,WidthRotorInlet];
            %
            %
            %% model gas stand combustion
            PressureTurbineOutlet = 101300.*ones(size(TurbineDataSort(:,1))); %Paa, assumed outlet pressure
            TemperatureTurbineInlet = (650+273).*ones(size(TurbineDataSort(:,1))); %K, assumed stagnation temperature at turbine inlet
            %combustion constants
            AirN2toO2 = 3.773;
            MW = [44.01,18.015,28.013,31.999]; %kg/kmol CO2,H2O,N2,O2
            RUniversal = 8314; %J/kmol/K
            RGas = RUniversal./MW; %J/kg/K CO2,H2O,N2,O2
            CpGasCoeffs = [1.0025e-7,-2.8008e-7,-1.9137e-7,-2.2590e-7;-6.2109e-4,7.9654e-4,4.9206e-4,4.4577e-4;1.2211,-0.0357,-0.1896,-0.0269;538.6045,1810.5,1055.8,901.3275]; %J/kg/K, rows are coeffs for T^3,T^2,T,constant [K], cols are CO2,H2O,N2,O2
            Air = [0,0,AirN2toO2,1]; %CO2,H2O,N2,O2
            MolFractionAir = Air./(Air*ones(4,1));
            MWAir = MolFractionAir*MW'; %kg/kmol
            MassFractionAir = MolFractionAir.*(MW./MWAir);
            %fuel
            Fuel = 2;
            if Fuel == 1 %diesel
                HeatingValueFuel = 43.2e6; %J/kg
                FuelHtoC = 1.8;
                MWFuel = 12+FuelHtoC;
                EquivalenceRatioEstimate = 25/14.6;
            elseif Fuel == 2 %natural gas
                HeatingValueFuel = 45e6; %J/kg
                FuelHtoC = 4;
                MWFuel = 12+FuelHtoC; %kg/kmol
                EquivalenceRatioEstimate = 2.3;
            end
            MolAirStoichiometric = (1+FuelHtoC/4).*Air; %kmol/kmol fuel
            MassAirStoichiometric = (MolAirStoichiometric.*MW)*ones(length(MW),1); %kg/kmol fuel
            FreshAirFuelRatioStoichiometric = MassAirStoichiometric/MWFuel;
            FreshAirFuelRatio = FreshAirFuelRatioStoichiometric.*EquivalenceRatioEstimate;
            %gas stand combustion
            MassFlowTurbinePhys = TurbineDataSort(:,3).*(PressureTurbineOutlet./1000)./TemperatureTurbineInlet.^0.5; %kg/sec
            MassFlowAir = MassFlowTurbinePhys./(1+1/FreshAirFuelRatio); %kg/sec
            MassFlowFuel = MassFlowAir./FreshAirFuelRatio; %kg/sec
            MolFlowFuel = MassFlowFuel./MWFuel; %kmol/sec
            MassFlowCombustionReactants = MassFlowAir*MassFractionAir; %kg/sec
            MolFlowCombustionReactants = MassFlowCombustionReactants./(ones(size(MassFlowAir))*MW); %kmol/sec
            MolFlowCombustionProducts = [MolFlowCombustionReactants(:,1)+MolFlowFuel,MolFlowCombustionReactants(:,2)+MolFlowFuel*(FuelHtoC/2),MolFlowCombustionReactants(:,3),MolFlowCombustionReactants(:,4)-MolFlowFuel*(1+FuelHtoC/4)]; %CO2,H2O,N2,O2
            MassFlowCombustionProducts = MolFlowCombustionProducts.*(ones(length(MassFlowAir),1)*MW);
            MassFractionExhaust = MassFlowCombustionProducts./(MassFlowTurbinePhys*ones(size(MW)));
            MolFractionExhaust = MolFlowCombustionProducts./((MolFlowCombustionProducts*ones(length(MW),1))*ones(size(MW)));
            CpTurbine = (([TemperatureTurbineInlet.^3,TemperatureTurbineInlet.^2,TemperatureTurbineInlet,ones(size(TemperatureTurbineInlet))]*CpGasCoeffs).*MassFractionExhaust)*ones(length(MW),1); %J/kg/K
            TemperatureTurbineInletModel = (MassFlowFuel.*HeatingValueFuel)./(MassFlowTurbinePhys.*CpTurbine)-273;
            RTurbine = MassFractionExhaust*RGas';
            CvTurbine = CpTurbine-RTurbine; %J/kg/K
            GammaTurbine = CpTurbine./CvTurbine;
            %solve additional inputs from mapping data
            SpeedTurbinePhys = TurbineDataSort(:,2).*TemperatureTurbineInlet.^0.5; %RPM
            SpeedRotorTip = SpeedTurbinePhys./60.*(2*pi*RadiusRotor); %m/sec
            EnthalpySpecificTurbineIdeal = CpTurbine.*TemperatureTurbineInlet.*(1-TurbineDataSort(:,4).^((1-GammaTurbine)./GammaTurbine));
            TemperatureTurbineOutletIdeal = TemperatureTurbineInlet-EnthalpySpecificTurbineIdeal./CpTurbine; %K
            CpTurbineOutletIdeal = (([TemperatureTurbineOutletIdeal.^3,TemperatureTurbineOutletIdeal.^2,TemperatureTurbineOutletIdeal,ones(size(TemperatureTurbineOutletIdeal))]*CpGasCoeffs).*MassFractionExhaust)*ones(length(MW),1); %J/kg/K
            CvTurbineOutletIdeal = CpTurbineOutletIdeal-RTurbine;
            GammaTurbineOutletIdeal = CpTurbineOutletIdeal./CvTurbineOutletIdeal;
            EnthalpySpecificTurbine = TurbineDataSort(:,5).*EnthalpySpecificTurbineIdeal;
            TemperatureTurbineOutlet = TemperatureTurbineInlet-EnthalpySpecificTurbine./CpTurbine; %K
            CpTurbineOutlet = (([TemperatureTurbineOutlet.^3,TemperatureTurbineOutlet.^2,TemperatureTurbineOutlet,ones(size(TemperatureTurbineOutlet))]*CpGasCoeffs).*MassFractionExhaust)*ones(length(MW),1); %J/kg/K
            CvTurbineOutlet = CpTurbineOutlet-RTurbine;
            GammaTurbineOutlet = CpTurbineOutlet./CvTurbineOutlet;
            UC0 = SpeedRotorTip./(2.*CpTurbine.*TemperatureTurbineInlet.*(1-TurbineDataSort(:,4).^((1-GammaTurbine)./GammaTurbine))).^0.5;
            TurbineDataSort = [TurbineDataSort(:,1:5),TemperatureTurbineInlet,PressureTurbineOutlet,SpeedRotorTip,CpTurbine,RTurbine,GammaTurbine,TemperatureTurbineOutletIdeal,...
                EnthalpySpecificTurbineIdeal,CpTurbineOutletIdeal,GammaTurbineOutletIdeal,TemperatureTurbineOutlet,EnthalpySpecificTurbine,CpTurbineOutlet,GammaTurbineOutlet,UC0];
            %
            %
            %% sort data and characterize data
            %find vane position breakpoints in input data
            VTGVanePosition = TurbineDataSort(:,1);
            MaxVTGVanePosition = max(VTGVanePosition);
            VTGVanePositionIndex = 1;
            for k = 2:length(VTGVanePosition)
                if abs(VTGVanePosition(k)-VTGVanePosition(k-1)) > 0.03*MaxVTGVanePosition
                    VTGVanePositionIndex = [VTGVanePositionIndex,k];
                end
            end
            VTGVanePositionIndex = [VTGVanePositionIndex,length(VTGVanePosition)+1];
            %find and sort speed breakpoints in vane position data
            for k1 = 1:length(VTGVanePositionIndex)-1
                TurbineDataVP = TurbineDataSort(VTGVanePositionIndex(k1):VTGVanePositionIndex(k1+1)-1,:);
                TurbineDataVPSpeedSort = sortrows(TurbineDataVP,2); %sort vane position data by corrected speed
                %find speedlines
                SpeedBreakpoints = 1;
                for k = 2:length(TurbineDataVPSpeedSort(:,2))
                    if abs(TurbineDataVPSpeedSort(k,2)-TurbineDataVPSpeedSort(k-1,2))/TurbineDataVPSpeedSort(k-1,2) > 0.05 %5% speed difference
                        SpeedBreakpoints = [SpeedBreakpoints,k];
                    end
                end
                SpeedBreakpoints = [SpeedBreakpoints,length(TurbineDataVPSpeedSort(:,1))+1];
                %sort speedline data by pressure ratio
                TurbineDataVPSort = zeros(size(TurbineDataVPSpeedSort));
                MaxEfficiency = zeros(length(SpeedBreakpoints)-1,3);
                for k = 1:length(SpeedBreakpoints)-1
                    TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,:) = sortrows(TurbineDataVPSpeedSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,:),4);
                    for k2 = SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1
                        if TurbineDataVPSort(k2,5) == max(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,5))
                            MaxEfficiencyIndex = k2;
                        end
                    end
                    MaxEfficiency(k,:) = [mean(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,2)),...
                        max(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,5)),MaxEfficiencyIndex];%+VTGVanePositionIndex(k1)-1];
                end
                %find choked flow data
                TurbineDataVPFlowSort = sortrows(TurbineDataVPSort,-3); %sort vane position data by corrected flow
                MassFlowCorrectedMax = (ones(1,5)*TurbineDataVPFlowSort(1:5,[3:4,6,10,11]))./5; %average flow,pressure ratio, inlet temp, R, gamma for 5 highest corrected flow points
                ChokeFlowMultiplier = [1,0.98,0.96,0.94,0.92,0.9];
                for k3 = 1:length(ChokeFlowMultiplier)
                    ChokeIndexVP = zeros(size(TurbineDataVPSort(:,1)));
                    for k = 1:length(TurbineDataVPSort(:,1))
                        if TurbineDataVPSort(k,4) > MassFlowCorrectedMax(2) && TurbineDataVPSort(k,3) < ChokeFlowMultiplier(k3)*MassFlowCorrectedMax(1)
                            ChokeIndexVP(k) = 1;
                        end
                    end
                    if k3 == 1
                        ChokeIndexGuess = ChokeIndexVP;
                    else %only add the next option if different choke points are identified
                        DeltaChoke = ones(1,length(ChokeIndexVP))*(ChokeIndexVP-ChokeIndexGuess(end,:));
                        if abs(DeltaChoke) > 0
                            ChokeIndexGuess = [ChokeIndexGuess,ChokeIndexVP];
                        end
                    end
                end
                %find heat transfer affected data points
                %assumed compressor operating condition
                TemperatureCompressorDischargeForHeatTransferError = 26+273; %K
                TemperatureCompressorInlet = 25+273; %K
                CpAir = 1005; %J/kg/K
                EnthalpyMaxForHeatTransferError = CpAir*(TemperatureCompressorDischargeForHeatTransferError-TemperatureCompressorInlet); %J/kg
                HeatTransferIndexVP = zeros(size(TurbineDataVPSort(:,1)));
                for k = 1:length(TurbineDataVPSort(:,17))
                    if TurbineDataVPSort(k,17) < EnthalpyMaxForHeatTransferError
                        HeatTransferIndexVP(k) = 1;
                    end
                end
                %assemble turbine map data for vane position
                TurbineData.VanePosition = mean(TurbineDataVPSort(:,1));
                TurbineData.Data = TurbineDataVPSort;
                TurbineData.MaxEfficiency = MaxEfficiency;
                TurbineData.ChokeIndexGuess = ChokeIndexGuess;
                TurbineData.HeatTransferIndex = HeatTransferIndexVP;
            end
            %
            %
            %% assemble output
            if isempty(obj.TurbineMapData)
                obj.TurbineMapData.TurbineMap = PathTurbineData;
                obj.TurbineMapData.DataType = DataType;
                obj.TurbineMapData.GeometryInput = GeometryTurbineData;
                obj.TurbineMapData.TurbineData = TurbineData;
            else
                obj.TurbineMapData(end+1).TurbineMap = PathTurbineData;
                obj.TurbineMapData(end).DataType = DataType;
                obj.TurbineMapData(end).GeometryInput = GeometryTurbineData;
                obj.TurbineMapData(end).TurbineData = TurbineData;
            end
            %
            %
        end %TurbineMapDataPrep
        
        function obj = BuildTurbineModel(obj)
            %fits turbine model to training data contained in
            %obj.TurbineMapData
            %
            global GeometryTurbine TurbineDataVPSortGoodData
            %
            %
            %% determine type of turbine
            NumberVanePositions = zeros(length(obj.TurbineMapData),1);
            NumberChokeAttempts = zeros(length(obj.TurbineMapData),1);
            TrainingIndex = [];
            for k = 1:length(obj.TurbineMapData)
                if strcmp(obj.TurbineMapData(k).DataType,'Training')
                    TrainingIndex = [TrainingIndex,k];
                end
                NumberVanePositions(k) = length(obj.TurbineMapData(k).TurbineData);
                NumberChokeAttempts(k) = length(obj.TurbineMapData(k).TurbineData(1).ChokeIndexGuess(1,:));
            end
            MeanVanePositions = (ones(1,length(TrainingIndex))*NumberVanePositions(TrainingIndex))/length(TrainingIndex);
            if MeanVanePositions == 1
                FixedGeometryFlag = 1;
            end
            %
            %
            if FixedGeometryFlag == 1 %make sure inputs are fixed goemetry before combining training maps
                %exhaustive search for choked flow - build model for each choked flow assumption
                NumberAttempts = max(NumberChokeAttempts(TrainingIndex));
                GeometryTurbineChokeSweep = zeros(NumberAttempts,6);
                CoeffsEnthalpyLossModelChokeSweep = zeros(NumberAttempts,3);
                CoeffsAngleIncidenceVsSpeedEfficiencyPeakChokeSweep = zeros(2*NumberAttempts,2);
                MachNoChokeSweep = zeros(NumberAttempts,1);
                SpeedRotorTipEfficiencyMaxChokeSweep = zeros(NumberAttempts,1);
                ErrorTrainingDataChokeSweep = zeros(2*NumberAttempts,9);
%                 TurbineDataModelChokeSweep = zeros(NumberAttempts,6*NumberAttempts);
                TurbineDataModelChokeSweep = [];
                %% apply geometry and combine training data sets
                GeometryTraining = [];
                GeometryTurbineModel = [];
                TrainingDataBreakpoints = 1;
                for k = 1:length(TrainingIndex)
                    GeometryTurbineModelTrainingData = [obj.TurbineMapData(TrainingIndex(k)).GeometryInput(1:3),...
                        obj.TurbineMapData(TrainingIndex(k)).GeometryInput(3)/obj.TurbineMapData(1).GeometryInput(3)];
                    GeometryTraining = [GeometryTraining;ones(size(obj.TurbineMapData(TrainingIndex(k)).TurbineData(1).Data(:,1)))*GeometryTurbineModelTrainingData];
                    GeometryTurbineModel = [GeometryTurbineModel;GeometryTurbineModelTrainingData];
                    TrainingDataBreakpoints = [TrainingDataBreakpoints,length(obj.TurbineMapData(TrainingIndex(k)).TurbineData(1).Data(:,1))+TrainingDataBreakpoints(end)];
                end
                %
                %attempt choked data point assumptions
                for k3 = 1:NumberAttempts
                    TrainingData = [];
                    for k1 = 1:length(TrainingIndex)
                        TrainingData = [TrainingData;obj.TurbineMapData(TrainingIndex(k1)).TurbineData(1).Data];
                        %find data points to use for model training
                        HeatTransferIndex = obj.TurbineMapData(TrainingIndex(k1)).TurbineData(1).HeatTransferIndex;
                        if length(obj.TurbineMapData(TrainingIndex(k1)).TurbineData(1).ChokeIndexGuess(1,:)) < k3
                            ChokeIndexGuess = obj.TurbineMapData(TrainingIndex(k1)).TurbineData(1).ChokeIndexGuess(:,end);
                        else
                            ChokeIndexGuess = obj.TurbineMapData(TrainingIndex(k1)).TurbineData(1).ChokeIndexGuess(:,k3);
                        end
                        %
                        %find data for fitting
                        HeatTransferChokeData = HeatTransferIndex+ChokeIndexGuess;
                        GoodDataIndex = [];
                        for k = 1:length(HeatTransferChokeData)
                            if HeatTransferChokeData(k) == 0
                                GoodDataIndex = [GoodDataIndex;k];
                            end
                        end
                    end
                    %
                    %data to use for fitting geometry and losses
                    GeometryTurbine = GeometryTraining(GoodDataIndex,:);
                    TurbineDataVPSortGoodData = TrainingData(GoodDataIndex,:);
                    %
                    %
                    %% find simplified geometry
                    RadiusRotorMin = min(GeometryTraining(:,1));
                    RadiusRotorMax = max(GeometryTraining(:,1));
                    RadiusRotorOutletMin = min(GeometryTraining(:,2));
                    RadiusRotorOutletMax = max(GeometryTraining(:,2));
                    AreaRadiusRatioVoluteThroat = obj.TurbineMapData(TrainingIndex(1)).GeometryInput(3); %use first turbine map for initial guesses
                    RadiusNoseRotorOutlet = obj.TurbineMapData(TrainingIndex(1)).GeometryInput(4); %use first turbine map for initial guesses
                    AngleBladeRotorOutlet = obj.TurbineMapData(TrainingIndex(1)).GeometryInput(5); %use first turbine map for initial guesses
                    %
                    GeometryLowerBounds = [0.1*AreaRadiusRatioVoluteThroat,0.1*RadiusRotorOutletMin,20/360*2*pi];
                    GeometryUpperBounds = [4*AreaRadiusRatioVoluteThroat,0.75*RadiusRotorOutletMax,80/360*2*pi];
                    DesignVariables = [AreaRadiusRatioVoluteThroat,RadiusNoseRotorOutlet,AngleBladeRotorOutlet];
                    GeometryTurbineFromFlow = fmincon(@FitTurbineGeometry7,DesignVariables,[],[],[],[],GeometryLowerBounds,GeometryUpperBounds);
                    GeometryTurbine(:,3) = GeometryTurbineFromFlow(1);%.*GeometryTurbine(:,4);
                    GeometryTurbine(:,5:6) = ones(size(GeometryTurbine(:,1)))*GeometryTurbineFromFlow(2:3);
                    %
                    %
                    %% find loss model
                    LossLowerBounds = [0.1*RadiusRotorMin,0,0,0];
                    LossUpperBounds = [0.4*RadiusRotorMax,20,10,10];
                    LossVariables = [0.25*RadiusRotorMin,1,0.3,0.1];
                    LossModelFromEfficiency = fmincon(@FitTurbineLoss7,LossVariables,[],[],[],[],LossLowerBounds,LossUpperBounds);
                    GeometryTurbine(:,7) = ones(size(GeometryTurbine(:,1)))*LossModelFromEfficiency(1);
                    %
                    %
                    %% fit choke model and optimum incidence model
                    [CoeffsAngleIncidenceVsSpeedEfficiencyPeak,MachNoChoke,SpeedRotorTipEfficiencyMax] = FitTurbineChoke7(GeometryTurbine,TurbineDataVPSortGoodData);
                    %
                    %
                    %% solve model for input data
                    GeometryTurbineFit = [GeometryTraining(:,1:2),ones(size(GeometryTraining(:,1))).*GeometryTurbine(1,3),GeometryTraining(:,4),ones(size(GeometryTraining(:,1)))*GeometryTurbine(1,5:7)];
                    TurbineDataInput = TrainingData(:,[8,4,6,9,10,11]);
                    TurbineDataModelAll = SolveTurbine7(GeometryTurbineFit,CoeffsAngleIncidenceVsSpeedEfficiencyPeak,LossModelFromEfficiency(2:end)',MachNoChoke,SpeedRotorTipEfficiencyMax,TurbineDataInput,obj.ConvergenceEfficiencyError);
                    TurbineDataModelAll(:,7) = TurbineDataInput(:,1)./(2.*TurbineDataInput(:,4).*TurbineDataInput(:,3).*(1-TurbineDataInput(:,2).^((1-TurbineDataInput(:,6))./TurbineDataInput(:,6)))).^0.5; %U/C0
                    ConvergenceIndexTraining = [];
                    for k = 1:length(TurbineDataModelAll(:,5))
                        if TurbineDataModelAll(k,5) <= obj.ConvergenceEfficiencyError
                            ConvergenceIndexTraining = [ConvergenceIndexTraining;k];
                        end
                    end
%                     TurbineDataModel = TurbineDataModelAll(ConvergenceIndexTraining,:);
                    ErrorTraining = CalculateError(TrainingData,TurbineDataModelAll,ConvergenceIndexTraining);
                    %
                    %
                    %% add choke guess to solution set
                    GeometryTurbineChokeSweep(k3,:) = GeometryTurbine(1,[1:3,5:7]);
                    CoeffsEnthalpyLossModelChokeSweep(k3,:) = LossModelFromEfficiency(2:end);
                    CoeffsAngleIncidenceVsSpeedEfficiencyPeakChokeSweep(2*k3-1:2*k3,:) = CoeffsAngleIncidenceVsSpeedEfficiencyPeak;
                    MachNoChokeSweep(k3) = MachNoChoke;
                    SpeedRotorTipEfficiencyMaxChokeSweep(k3) = SpeedRotorTipEfficiencyMax;
                    ErrorTrainingDataChokeSweep(2*k3-1:2*k3,:) = ErrorTraining;
%                     TurbineDataModelChokeSweep(:,(6*(k3-1)+1):(6*k3)) = TurbineDataModelAll;
%                     TurbineDataModelChokeSweep = [TurbineDataModelChokeSweep,TurbineDataModelAll];
                end %number of choke attempts
                %
                %
                %% pick best fit choke index based on efficiency RMS error
                ErrorChokeSweep = zeros(1,NumberAttempts);
                for k3 = 1:NumberAttempts
                    ErrorChokeSweep(k3) = ErrorTrainingDataChokeSweep(2*k3,2);
                end
                MinErrorChokeSweep = min(ErrorChokeSweep);
                for k3 = 1:NumberAttempts
                    if ErrorChokeSweep(k3) == MinErrorChokeSweep
                        GeometryTurbineFinal = GeometryTurbineChokeSweep(k3,:);
                        CoeffsEnthalpyLossModel = CoeffsEnthalpyLossModelChokeSweep(k3,:);
                        CoeffsAngleIncidenceVsSpeedEfficiencyPeak = CoeffsAngleIncidenceVsSpeedEfficiencyPeakChokeSweep(2*k3-1:2*k3,:);
                        MachNoChoke = MachNoChokeSweep(k3);
                        SpeedRotorTipEfficiencyMax = SpeedRotorTipEfficiencyMaxChokeSweep(k3);
                        ErrorTraining = ErrorTrainingDataChokeSweep(2*k3-1:2*k3,:);
%                         TurbineDataModel = TurbineDataModelChokeSweep(:,(7*(k3-1)+1):(7*k3));
                    end
                end
                %
                %
                %% write results to object
                GeometryTurbineUpdate = [GeometryTurbineFinal(1:2),obj.TurbineMapData(TrainingIndex(1)).GeometryInput(3),GeometryTurbineFinal(3:end)];
                TurbineModelTable = array2table([GeometryTurbineUpdate,CoeffsEnthalpyLossModel,CoeffsAngleIncidenceVsSpeedEfficiencyPeak(:,1)',MachNoChoke]);
                TurbineModelTable.Properties.VariableNames = {'RotorRadius_m','RotorOutletRadius_m','AR_Geometry_m','AR_Learned_m','RotorOutletNoseRadius_m','BladeOutletAngle_rad','RotorInletWidth_m',...
                    'IncidenceLossCoefficient','RotorFrictionCoefficient','StatorFrictionCoefficient','OptimumIncidenceVsSpeed_1','OptimumIncidenceVsSpeed_0','MachNoChokeModel'};
                TurbineScalingModelTable = array2table([CoeffsAngleIncidenceVsSpeedEfficiencyPeak(:,2)',SpeedRotorTipEfficiencyMax]);
                TurbineScalingModelTable.Properties.VariableNames = {'OptimumFlowCoefficientVsSpeed_1','OptimumFlowCoefficientVsSpeed_0','RotorSpeedForMaxEfficiency'};
                ErrorTrainingTable = array2table(ErrorTraining);
                ErrorTrainingTable.Properties.VariableNames = {'MeanError','ErrorStDev','RMSError','MeanError_Unchoke','ErrorStDev_Unchoke','RMSError_Unchoke','MeanError_LowFlow','ErrorStDev_LowFlow','RMSError_LowFlow'};
                ErrorTrainingTable.Properties.RowNames = {'Flow','Efficiency'};
                %
                obj.TurbineModel.Model = TurbineModelTable;
                obj.TurbineModel.Scale = TurbineScalingModelTable;
                obj.TurbineModel.ErrorTrainingData = ErrorTrainingTable;
%                 obj.TurbineModel.TurbineMappingModel = TurbineMappingModel;
%
            end %if fixed geometry
        end %BuildTurbineModel
        
        function obj = SolveTurbineModel(obj)
            %
            %
            %extract turbine model parameters
            TurbineModelInput = table2array(obj.TurbineModel.Model);
            TurbineModelScale = table2array(obj.TurbineModel.Scale);
            GeometryTurbine = [TurbineModelInput([1,2,4]),1,TurbineModelInput(5:7)];
            CoeffsAngleIncidenceVsSpeedEfficiencyPeak = [TurbineModelInput(11:12)',TurbineModelScale(1:2)'];
            CoeffsEnthalpyLossModel = TurbineModelInput(8:10);
            MachNoChoke = TurbineModelInput(13);
            SpeedRotorTipEfficiencyMax = TurbineModelScale(3);
            %
            %
            %find first training data set
            TrainingIndex = [];
            for k = 1:length(obj.TurbineMapData)
                if strcmp(obj.TurbineMapData(k).DataType,'Training')
                    TrainingIndex = [TrainingIndex,k];
                end
            end
            %
            for k1 = 1:length(obj.TurbineMapData)
                TurbineDataVPSort = obj.TurbineMapData(k1).TurbineData(1).Data;
                ScaledMapFlag = 0;
                if TurbineDataVPSort == 0
                    ScaledMapFlag = 1;
                    TurbineDataVPSort = obj.TurbineMapData(TrainingIndex(1)).TurbineData(1).Data;
                end
                obj.TurbineModel.TurbineMappingModel(k1).TurbineMap = obj.TurbineMapData(k1).TurbineMap;
                obj.TurbineModel.TurbineMappingModel(k1).DataType = obj.TurbineMapData(k1).DataType;
                obj.TurbineModel.TurbineMappingModel(k1).ScaleFactors = [obj.TurbineMapData(k1).GeometryInput(1:2)./obj.TurbineMapData(TrainingIndex(1)).GeometryInput(1:2),1,...
                    obj.TurbineMapData(k1).GeometryInput(3)/obj.TurbineMapData(TrainingIndex(1)).GeometryInput(3),obj.TurbineMapData(k1).GeometryInput(1)/obj.TurbineMapData(TrainingIndex(1)).GeometryInput(1),...
                    1,obj.TurbineMapData(k1).GeometryInput(1)/obj.TurbineMapData(TrainingIndex(1)).GeometryInput(1)];
                GeometryTurbineScale = GeometryTurbine.*obj.TurbineModel.TurbineMappingModel(k1).ScaleFactors;
                %
                %
                %% find speedlines in raw data
                SpeedBreakpoints = 1;
                for k = 2:length(TurbineDataVPSort(:,2))
                    if abs(TurbineDataVPSort(k,2)-TurbineDataVPSort(k-1,2))/TurbineDataVPSort(k-1,2) > 0.05 %5% speed difference
                        SpeedBreakpoints = [SpeedBreakpoints,k];
                    end
                end
                SpeedBreakpoints = [SpeedBreakpoints,length(TurbineDataVPSort(:,1))+1];
                %
                %
                %% solve turbine map data points
                if ScaledMapFlag == 0
                    PressureRatio = TurbineDataVPSort(:,4);
                    SpeedRotorTip = TurbineDataVPSort(:,8);
                    TemperatureTurbineInlet = TurbineDataVPSort(:,6);
                    CpTurbine = TurbineDataVPSort(:,9);
                    RTurbine = TurbineDataVPSort(:,10);
                    GammaTurbine = TurbineDataVPSort(:,11);
                    TurbineDataInput = [SpeedRotorTip,PressureRatio,TemperatureTurbineInlet,CpTurbine,RTurbine,GammaTurbine];
                    GeometryTurbineScaleFinal = ones(size(SpeedRotorTip))*GeometryTurbineScale;
                    TurbineDataModelAll = SolveTurbine7(GeometryTurbineScaleFinal,CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsEnthalpyLossModel',MachNoChoke,SpeedRotorTipEfficiencyMax,TurbineDataInput,obj.ConvergenceEfficiencyError);
                    UC0ModelAll = TurbineDataInput(:,1)./(2.*TurbineDataInput(:,4).*TurbineDataInput(:,3).*(1-TurbineDataInput(:,2).^((1-TurbineDataInput(:,6))./TurbineDataInput(:,6)))).^0.5;
                    ConvergenceIndex = [];
                    for k = 1:length(TurbineDataModelAll(:,5))
                        if TurbineDataModelAll(k,5) <= 1.5*obj.ConvergenceEfficiencyError% && HeatTransferIndex(k) == 0
                            ConvergenceIndex = [ConvergenceIndex;k];
                        end
                    end
                    TurbineDataModel = TurbineDataModelAll(ConvergenceIndex,:);
                    UC0Model = UC0ModelAll(ConvergenceIndex,:);
                    obj.TurbineModel.TurbineMappingModel(k1).ModelOutput = [obj.TurbineMapData(k1).TurbineData(1).VanePosition.*ones(size(TurbineDataModel(:,1))),TurbineDataModel,UC0Model];
                    ErrorTraining = CalculateError(TurbineDataVPSort,TurbineDataModelAll,ConvergenceIndex);
                    ErrorTrainingTable = array2table(ErrorTraining);
                    ErrorTrainingTable.Properties.VariableNames = {'MeanError','ErrorStDev','RMSError','MeanError_Unchoke','ErrorStDev_Unchoke','RMSError_Unchoke','MeanError_LowFlow','ErrorStDev_LowFlow','RMSError_LowFlow'};
                    ErrorTrainingTable.Properties.RowNames = {'Flow','Efficiency'};
                    obj.TurbineModel.TurbineMappingModel(k1).ModelError = ErrorTrainingTable;
                end
                %
                %
                %% extend speedlines for training data model and solve
                PressureRatioModel = [];
                SpeedRotorTipModel = [];
                TemperatureTurbineInletModel = [];
                CpTurbineModel = [];
                RTurbineModel = [];
                GammaTurbineModel = [];
                PressureRatioSpeedline = (1:0.02:5)';
                %extend speedlines in mapping data
                for k = 1:length(SpeedBreakpoints)-1
                    SpeedRotorTipSpeedline = (mean(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,8))).*ones(size(PressureRatioSpeedline));
                    TemperatureTurbineInletSpeedline = (mean(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,6))).*ones(size(PressureRatioSpeedline));
                    CpTurbineSpeedline = (mean(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,9))).*ones(size(PressureRatioSpeedline));
                    RTurbineSpeedline = (mean(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,10))).*ones(size(PressureRatioSpeedline));
                    GammaTurbineSpeedline = (mean(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,11))).*ones(size(PressureRatioSpeedline));
                    PressureRatioModel = [PressureRatioModel;PressureRatioSpeedline];
                    SpeedRotorTipModel = [SpeedRotorTipModel;SpeedRotorTipSpeedline];
                    TemperatureTurbineInletModel = [TemperatureTurbineInletModel;TemperatureTurbineInletSpeedline];
                    CpTurbineModel = [CpTurbineModel;CpTurbineSpeedline];
                    RTurbineModel = [RTurbineModel;RTurbineSpeedline];
                    GammaTurbineModel = [GammaTurbineModel;GammaTurbineSpeedline];
                end
                %add high speedlines
                SpeedRotorTipMaxModel = 500; %m/sec, update with user input
                MaxSpeedRotorTip = max(SpeedRotorTipModel);
                if MaxSpeedRotorTip < SpeedRotorTipMaxModel-50
                    SpeedRotorTipHighSpeedModel = ((linspace(MaxSpeedRotorTip^2,SpeedRotorTipMaxModel^2,3)).^0.5);
                    for k = 2:length(SpeedRotorTipHighSpeedModel)
                        SpeedRotorTipSpeedline = SpeedRotorTipHighSpeedModel(k).*ones(size(PressureRatioSpeedline));
                        TemperatureTurbineInletSpeedline = TemperatureTurbineInletModel(end).*ones(size(PressureRatioSpeedline));
                        CpTurbineSpeedline = CpTurbineModel(end).*ones(size(PressureRatioSpeedline));
                        RTurbineSpeedline = RTurbineModel(end).*ones(size(PressureRatioSpeedline));
                        GammaTurbineSpeedline = GammaTurbineModel(end).*ones(size(PressureRatioSpeedline));
                        PressureRatioModel = [PressureRatioModel;PressureRatioSpeedline];
                        SpeedRotorTipModel = [SpeedRotorTipModel;SpeedRotorTipSpeedline];
                        TemperatureTurbineInletModel = [TemperatureTurbineInletModel;TemperatureTurbineInletSpeedline];
                        CpTurbineModel = [CpTurbineModel;CpTurbineSpeedline];
                        RTurbineModel = [RTurbineModel;RTurbineSpeedline];
                        GammaTurbineModel = [GammaTurbineModel;GammaTurbineSpeedline];
                    end
                end                
                %add low speedlines
                SpeedRotorTipMinModel = 100; %m/sec, update with user input
                MinSpeedRotorTip = min(SpeedRotorTipModel);
                if MinSpeedRotorTip > SpeedRotorTipMinModel+50
                    SpeedRotorTipLowSpeedModel = ((linspace(SpeedRotorTipMinModel^2,MinSpeedRotorTip^2,3)).^0.5);
                    for k = 1:length(SpeedRotorTipLowSpeedModel)-1
                        SpeedRotorTipSpeedline = SpeedRotorTipLowSpeedModel(k).*ones(size(PressureRatioSpeedline));
                        TemperatureTurbineInletSpeedline = TemperatureTurbineInletModel(1).*ones(size(PressureRatioSpeedline));
                        CpTurbineSpeedline = CpTurbineModel(1).*ones(size(PressureRatioSpeedline));
                        RTurbineSpeedline = RTurbineModel(1).*ones(size(PressureRatioSpeedline));
                        GammaTurbineSpeedline = GammaTurbineModel(1).*ones(size(PressureRatioSpeedline));
                        PressureRatioModel = [PressureRatioModel;PressureRatioSpeedline];
                        SpeedRotorTipModel = [SpeedRotorTipModel;SpeedRotorTipSpeedline];
                        TemperatureTurbineInletModel = [TemperatureTurbineInletModel;TemperatureTurbineInletSpeedline];
                        CpTurbineModel = [CpTurbineModel;CpTurbineSpeedline];
                        RTurbineModel = [RTurbineModel;RTurbineSpeedline];
                        GammaTurbineModel = [GammaTurbineModel;GammaTurbineSpeedline];
                    end
                end
                %solve extended map
                TurbineDataInput = [SpeedRotorTipModel,PressureRatioModel,TemperatureTurbineInletModel,CpTurbineModel,RTurbineModel,GammaTurbineModel];
                GeometryTurbineScaleFinal = ones(size(SpeedRotorTipModel))*GeometryTurbineScale;
                TurbineDataModelAll = SolveTurbine7(GeometryTurbineScaleFinal,CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsEnthalpyLossModel',MachNoChoke,SpeedRotorTipEfficiencyMax,TurbineDataInput,obj.ConvergenceEfficiencyError);
                UC0ModelAll = TurbineDataInput(:,1)./(2.*TurbineDataInput(:,4).*TurbineDataInput(:,3).*(1-TurbineDataInput(:,2).^((1-TurbineDataInput(:,6))./TurbineDataInput(:,6)))).^0.5;
                ConvergenceIndex = [];
                for k = 1:length(TurbineDataModelAll(:,5))
                    if TurbineDataModelAll(k,5) <= 1.5*obj.ConvergenceEfficiencyError% && HeatTransferIndex(k) == 0
                        ConvergenceIndex = [ConvergenceIndex;k];
                    end
                end
                TurbineDataModel = TurbineDataModelAll(ConvergenceIndex,:);
                UC0Model = UC0ModelAll(ConvergenceIndex,:);
                obj.TurbineModel.TurbineMappingModel(k1).ModelExtended = [obj.TurbineMapData(k1).TurbineData(1).VanePosition.*ones(size(TurbineDataModel(:,1))),TurbineDataModel,UC0Model];
                %
                %
            end
        end %SolveTurbineModel
        
        function obj = CreateScaledTurbine(obj,ScaleGeometry)
            %
            %
            %find first training data set
            TrainingIndex = [];
            for k = 1:length(obj.TurbineMapData)
                if strcmp(obj.TurbineMapData(k).DataType,'Training')
                    TrainingIndex = [TrainingIndex,k];
                end
            end
            %
            ScaleData.VanePosition = 1;
            ScaleData.Data = 0;
            ScaleData.MaxEfficiency = 0;
            ScaleData.ChokeIndexGuess = 0;
            ScaleData.HeatTransferIndex = 0;
            %
            GeometryTraining = obj.TurbineMapData(TrainingIndex(1)).GeometryInput;
            GeometryTurbineScale = GeometryTraining;
            GeometryTurbineScale(1:3) = ScaleGeometry.*GeometryTraining(1:3);%
            obj.TurbineMapData(end+1).TurbineMap = 'NewScale';
            obj.TurbineMapData(end).DataType = 'ScaledMap';
            obj.TurbineMapData(end).GeometryInput = GeometryTurbineScale;
            obj.TurbineMapData(end).TurbineData = ScaleData;
            %
            %
        end
        
        function PlotModel(obj)
            %
            %
            ColorOrder = [0,0,0;1,0,0;0,1,0;0.85,0.35,0;1,0,1;0,1,1;0.25,0.25,0;0.25,0,0.25;0,0.25,0.25;0.5,0.5,0;0.5,0,0.5;0,0.5,0.5;0.75,0.75,0;0.75,0,0.75;0,0.75,0.75];
            ColorOrderLow = [0,0.375,0.375;0.375,0.375,0;0.375,0,0.375;0,0.625,0.625;0.625,0.625,0;0.625,0,0.625];
            ColorOrderHigh = [0,0.8,0.8;0.8,0.8,0;0.8,0,0.8;0,0.2,0.2;0.2,0.2,0;0.2,0,0.2];
            MarkerOrder = {'x','s','^','d','v','*','+','p','<','h','>'};
            %
            %
            for k1 = 1:length(obj.TurbineMapData)
                %initialize maps
                figure(10+k1)
                hold on
                xlabel('Pressure Ratio')
                ylabel('Corrected Mass Flow (kg/sec*K^1^/^2/kPaa)')
            %     axis([1,MaxPressureRatioTurbineFit,0,1.1*max(TurbineDataSort(:,3))])
                figure(20+k1)
                hold on
                xlabel('Pressure Ratio')
                ylabel('Efficiency')
                axis([1,5,0.3,0.9])
                figure(30+k1)
                hold on
                xlabel('U/C0')
                ylabel('Efficiency')
                axis([0,1.5,0.3,0.9])
                %% plot map data
                if strcmp(obj.TurbineMapData(k1).DataType,'Training') || strcmp(obj.TurbineMapData(k1).DataType,'Validation')
                    TurbineDataVPSort = obj.TurbineMapData(k1).TurbineData(1).Data;
                    SpeedBreakpoints = 1;
                    for k = 2:length(TurbineDataVPSort(:,2))
                        if abs(TurbineDataVPSort(k,2)-TurbineDataVPSort(k-1,2))/TurbineDataVPSort(k-1,2) > 0.05 %5% speed difference
                            SpeedBreakpoints = [SpeedBreakpoints,k];
                        end
                    end
                    SpeedBreakpoints = [SpeedBreakpoints,length(TurbineDataVPSort(:,1))+1];
                    %
                    %plot map data
                    for k = 1:length(SpeedBreakpoints)-1
    %                     SpeedLineLegend(k) = {[num2str(SpeedTipMean),' m/sec']};
                        figure(10+k1)
                        plot(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,4),TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,3),'-','Marker',char(MarkerOrder(k)),'Color',ColorOrder(k,:))
                        figure(20+k1)
                        plot(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,4),TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,5),'-','Marker',char(MarkerOrder(k)),'Color',ColorOrder(k,:))
                        figure(30+k1)
                        plot(TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,20),TurbineDataVPSort(SpeedBreakpoints(k):SpeedBreakpoints(k+1)-1,5),'-','Marker',char(MarkerOrder(k)),'Color',ColorOrder(k,:)) 
                    end
                end
                %
                %
                %% plot model
                TurbineDataModel = obj.TurbineModel.TurbineMappingModel(k1).ModelExtended;
                SpeedBreakpointsTrainingModel = 1;
                for k = 2:length(TurbineDataModel(:,2))
                    if abs(TurbineDataModel(k,2)-TurbineDataModel(k-1,2))/TurbineDataModel(k-1,2) > 0.05 %5% speed difference
                        SpeedBreakpointsTrainingModel = [SpeedBreakpointsTrainingModel,k];
                    end
                end
                SpeedBreakpointsTrainingModel = [SpeedBreakpointsTrainingModel,length(TurbineDataModel(:,1))+1];
                %plot model
                for k = 1:length(SpeedBreakpointsTrainingModel)-1
%                     SpeedLineLegend(end+1) = {[num2str(SpeedTipMean),' m/sec']};
                    figure(10+k1)
                    plot(TurbineDataModel(SpeedBreakpointsTrainingModel(k):SpeedBreakpointsTrainingModel(k+1)-1,4),TurbineDataModel(SpeedBreakpointsTrainingModel(k):SpeedBreakpointsTrainingModel(k+1)-1,3),'--','Color',ColorOrder(k,:))
                    figure(20+k1)
                    plot(TurbineDataModel(SpeedBreakpointsTrainingModel(k):SpeedBreakpointsTrainingModel(k+1)-1,4),TurbineDataModel(SpeedBreakpointsTrainingModel(k):SpeedBreakpointsTrainingModel(k+1)-1,5),'--','Color',ColorOrder(k,:))
                    figure(30+k1)
                    plot(TurbineDataModel(SpeedBreakpointsTrainingModel(k):SpeedBreakpointsTrainingModel(k+1)-1,8),TurbineDataModel(SpeedBreakpointsTrainingModel(k):SpeedBreakpointsTrainingModel(k+1)-1,5),'--','Color',ColorOrder(k,:))
                end
                
            end
        end %PlotModel
        
    end
end

