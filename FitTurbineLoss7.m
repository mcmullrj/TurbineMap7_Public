function Error2Efficiency = FitTurbineLoss7(LossModelVariables)
%loss model optimization function
%created by Bob McMullen
%
%
%% design variables
WidthRotorInlet = LossModelVariables(1);
CoeffsEnthalpyLossModel = LossModelVariables(2:end)';
%
%
%% turbine map and geometry inputs
global GeometryTurbine TurbineDataVPSortGoodData
% extract geometry
AreaRadiusRatioVoluteThroat = GeometryTurbine(:,3); %m
AngleBladeRotorOutlet = GeometryTurbine(:,6); %rad
RadiusRotor = GeometryTurbine(:,1); %m
RadiusTipRotorOutlet = GeometryTurbine(:,2); %m
RadiusNoseRotorOutlet = GeometryTurbine(:,5); %m
AreaRadiusRatioMultiplier = GeometryTurbine(:,4);
%extract turbine performance data
SpeedRotorTip = TurbineDataVPSortGoodData(:,8); %m/sec
MassFlowCorrected =  TurbineDataVPSortGoodData(:,3)./1000; %kg/sec*sqrtK/Paa
PressureRatioTurbine =  TurbineDataVPSortGoodData(:,4);
EfficiencyTurbine =  TurbineDataVPSortGoodData(:,5);
TemperatureTurbineInlet =  TurbineDataVPSortGoodData(:,6); %K
% PressureTurbineOutlet =  TurbineDataVPSortGoodData(:,7); %Paa
RTurbine =  TurbineDataVPSortGoodData(:,10); %J/kg/K
CpTurbine =  TurbineDataVPSortGoodData(:,9); %J/kg/K
GammaTurbine =  TurbineDataVPSortGoodData(:,11);
EnthalpySpecificTurbineIdeal =  TurbineDataVPSortGoodData(:,13); %J/kg
% TemperatureTurbineOutletIdeal =  TurbineDataVPSortGoodData(:,12); %K
% CpTurbineOutletIdeal =  TurbineDataVPSortGoodData(:,14); %J/kg/K
% GammaTurbineOutletIdeal =  TurbineDataVPSortGoodData(:,15); %J/kg/K
EnthalpySpecificTurbine =  TurbineDataVPSortGoodData(:,17); %J/kg
% TemperatureTurbineOutlet =  TurbineDataVPSortGoodData(:,16); %K
% CpTurbineOutlet =  TurbineDataVPSortGoodData(:,18); %J/kg/K
GammaTurbineOutlet =  TurbineDataVPSortGoodData(:,19); %J/kg/K
%
%
%% update turbine geometry
% RadiusRMSRotorOutlet = ((RadiusTipRotorOutlet^2+RadiusNoseRotorOutlet^2)/2)^0.5; %m
AreaRotorOutlet = pi.*(RadiusTipRotorOutlet.^2-RadiusNoseRotorOutlet.^2); %m^2
%
%
%% velocity triangles
% section 1 velocity triangle
%assumptions: isentropic flow from 0 to 1 for Vm1
MassFlowParameterRotor = MassFlowCorrected.*RTurbine.*TemperatureTurbineInlet.^0.5;
VelocityTangentialRotorInlet = MassFlowParameterRotor./(AreaRadiusRatioVoluteThroat.*AreaRadiusRatioMultiplier)./RadiusRotor; %m/sec
VelocityTangentialRelativeRotorInlet = VelocityTangentialRotorInlet-SpeedRotorTip;
%find M1
MachNoRotorInletGuess = (0.01:0.01:1);
VelocityTangentialRotorInlet2FromTriangle = ((GammaTurbine.*RTurbine.*TemperatureTurbineInlet)*MachNoRotorInletGuess.^2)./(1+((GammaTurbine-1)./2)*MachNoRotorInletGuess.^2)-...
    ((MassFlowParameterRotor./(2.*pi.*RadiusRotor.*WidthRotorInlet))*ones(size(MachNoRotorInletGuess))).*(1+((GammaTurbine-1)./2)*MachNoRotorInletGuess.^2).^((1./(GammaTurbine-1))*ones(size(MachNoRotorInletGuess)));
Error2VelocityTangentialRotorInlet2 = (VelocityTangentialRotorInlet2FromTriangle-VelocityTangentialRotorInlet.^2*ones(size(MachNoRotorInletGuess))).^2;
MachNoRotorInlet = zeros(size(VelocityTangentialRotorInlet));
for k = 1:length(VelocityTangentialRotorInlet)
    Error2Min = min(Error2VelocityTangentialRotorInlet2(k,:));
    for k1 = 1:length(MachNoRotorInletGuess)
        if Error2VelocityTangentialRotorInlet2(k,k1) == Error2Min
            MachNoRotorInlet(k) = MachNoRotorInletGuess(k1);
        end
    end
end
%solve for Vm1
VelocityMeridionalRotorInlet = MassFlowParameterRotor./(2.*pi.*RadiusRotor.*WidthRotorInlet).*(1+(GammaTurbine-1)./2.*MachNoRotorInlet.^2).^(1./(GammaTurbine-1)); %m/sec
VelocityRotorInlet = (VelocityTangentialRotorInlet.^2+VelocityMeridionalRotorInlet.^2).^0.5;
VelocityRelativeRotorInlet = (VelocityTangentialRelativeRotorInlet.^2+VelocityMeridionalRotorInlet.^2).^0.5; %m/sec
AngleRelativeRotorInlet = atan(VelocityTangentialRelativeRotorInlet./VelocityMeridionalRotorInlet); %rad
%
%section 2 velocity triangle
TemperatureRatioTurbine = 1-EnthalpySpecificTurbine./CpTurbine./TemperatureTurbineInlet;
VelocityMeridionalRotorOutlet = (MassFlowParameterRotor./AreaRotorOutlet./TemperatureRatioTurbine.*PressureRatioTurbine); %m/sec
% VelocityTangentialRotorOutlet = SpeedRotorTip.*(RadiusRMSRotorOutlet/RadiusRotor)-VelocityMeridionalRotorOutlet.*(1./cos(AngleBladeRotorOutlet).^2-1).^0.5; %m/sec
VelocityRelativeRotorOutlet = VelocityMeridionalRotorOutlet./cos(AngleBladeRotorOutlet); %m/sec
% %modify choked flow points
% TemperatureTurbineOutlet = TemperatureTurbineInlet./TemperatureRatioTurbine; %K
% MachNoRelativeRotorOutlet = VelocityRelativeRotorOutlet./(GammaTurbineOutlet.*RTurbine.*TemperatureTurbineOutlet).^0.5;
% MachNoChoke = max(MachNoRelativeRotorOutlet(GoodDataIndex));
% for k = 1:length(MachNoRelativeRotorOutlet)
%     if MachNoRelativeRotorOutlet(k) > MachNoChoke
%         TemperatureRelativeRotorOutlet = TemperatureTurbineOutlet(k)-(VelocityRelativeRotorOutlet(k)^2/2/CpTurbineOutlet(k)); %K
%         VelocityRelativeRotorOutlet(k) = (GammaTurbineOutlet(k)*RTurbine(k)*TemperatureRelativeRotorOutlet)^0.5; %m/sec
%         AngleFlowRotorOutlet = VelocityRelativeRotorOutlet
%     end
% end
%
%
%% fit peak efficiency incidence angle vs. tip speed
%find data from different training maps
AngleIncidenceEfficiencyPeak = [];
MapBreakpoints = 1;
for k = 2:length(AreaRadiusRatioMultiplier)
    if abs(AreaRadiusRatioMultiplier(k)-AreaRadiusRatioMultiplier(k-1)) > 0.04 %4% change in A/R ratio
        MapBreakpoints = [MapBreakpoints,k];
    end
end
MapBreakpoints = [MapBreakpoints,length(AreaRadiusRatioMultiplier)+1];
%find speedlines and fit peak efficiency for each data set
for k1 = 1:length(MapBreakpoints)-1
    %find speed breakpoints
    SpeedRotorTurbineMap = SpeedRotorTip(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
    EfficiencyTurbineMap = EfficiencyTurbine(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
    AngleRelativeRotorInletMap = AngleRelativeRotorInlet(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
    SpeedBreakpoints = 1;
    for k2 = 2:length(SpeedRotorTurbineMap)
        if abs(SpeedRotorTurbineMap(k2)-SpeedRotorTurbineMap(k2-1)) > 15 %3% of 500 m/sec
            SpeedBreakpoints = [SpeedBreakpoints,k2];
        end
    end
    SpeedBreakpoints = [SpeedBreakpoints,length(SpeedRotorTurbineMap)+1];
    %find peak efficiency
    MaxEfficiencyIndex = [];
    for k2 = 1:length(SpeedBreakpoints)-1
        MaxEfficiencySpeedline = max(EfficiencyTurbineMap(SpeedBreakpoints(k2):SpeedBreakpoints(k2+1)-1));
        for k3 = SpeedBreakpoints(k2):SpeedBreakpoints(k2+1)-1
            if EfficiencyTurbineMap(k3) == MaxEfficiencySpeedline
                MaxEfficiencyIndex = [MaxEfficiencyIndex;k3];
            end
        end
    end
    CoeffsAngleIncidenceVsSpeedEfficiencyPeak = [SpeedRotorTurbineMap(MaxEfficiencyIndex),ones(length(MaxEfficiencyIndex),1)]\AngleRelativeRotorInletMap(MaxEfficiencyIndex);
    AngleIncidenceEfficiencyPeak = [AngleIncidenceEfficiencyPeak;[SpeedRotorTurbineMap,ones(size(SpeedRotorTurbineMap))]*CoeffsAngleIncidenceVsSpeedEfficiencyPeak];
end
%
%
%% enthalpy loss models
%enthalpy loss due to incidence
EnthalpyLossModelIncidence = (sin(AngleRelativeRotorInlet-AngleIncidenceEfficiencyPeak).^2.*VelocityRelativeRotorInlet.^2./2)./EnthalpySpecificTurbineIdeal; %J/kg
%
%enthalpy loss due to flow friction
EnthalpyLossModelRotorFriction = (cos(AngleRelativeRotorInlet-AngleIncidenceEfficiencyPeak).^2.*VelocityRelativeRotorInlet.^2./2+VelocityRelativeRotorOutlet.^2./2)./EnthalpySpecificTurbineIdeal; %J/kg
%
%enthalpy loss due to bearing loss
% EnthalpyLossModelBearings = SpeedRotorTip.^2./(1e5.*MassFlowCorrected.*(PressureRatioTurbine.*101300)./TemperatureTurbineInlet.^0.5)./EnthalpySpecificTurbineIdeal;
%
%enthalpy loss stator
EnthalpyLossModelStator = VelocityRotorInlet.^2./2./EnthalpySpecificTurbineIdeal;
%
%overall loss model
% EnthalpyLossModel = [EnthalpyLossModelBearings,EnthalpyLossModelIncidence,EnthalpyLossModelRotorFriction,EnthalpyLossModelStator]*CoeffsEnthalpyLossModel;
EnthalpyLossModel = [EnthalpyLossModelIncidence,EnthalpyLossModelRotorFriction,EnthalpyLossModelStator]*CoeffsEnthalpyLossModel;
EfficiencyModel = 1-EnthalpyLossModel;
%
%
%% calculate efficiency model error
MaxEfficiencyTurbine = max(EfficiencyTurbine);
ErrorEfficiency = (EfficiencyModel-EfficiencyTurbine)./MaxEfficiencyTurbine; %normalize error
% figure(3)
% plot(PressureRatioTurbine,ErrorEfficiency,'o')
Error2Efficiency = ((ones(1,length(ErrorEfficiency))*ErrorEfficiency.^2)/length(ErrorEfficiency));
%
%
%
%
end