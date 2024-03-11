function [CoeffsAngleIncidenceVsSpeedEfficiencyPeak,MachNoChoke,SpeedRotorTipEfficiencyMax] = FitTurbineChoke7(GeometryTurbine,TurbineDataVPSortGoodData)
%fit additional model inputs for optimized geometry and loss models
%created by Bob McMullen
%
%
%% turbine map and geometry inputs
% extract geometry
AreaRadiusRatioVoluteThroat = GeometryTurbine(:,3); %m
AreaRadiusRatioMultiplier = GeometryTurbine(:,4);
AngleBladeRotorOutlet = GeometryTurbine(:,6); %rad
RadiusRotor = GeometryTurbine(:,1); %m
RadiusTipRotorOutlet = GeometryTurbine(:,2); %m
RadiusNoseRotorOutlet = GeometryTurbine(:,5); %m
WidthRotorInlet = GeometryTurbine(:,7); %m
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
% EnthalpySpecificTurbineIdeal =  TurbineDataVPSortGoodData(:,13); %J/kg
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
AreaRotorOutlet = pi.*(RadiusTipRotorOutlet.^2-RadiusNoseRotorOutlet.^2); %m^2
% RadiusRMSRotorOutlet = ((RadiusTipRotorOutlet^2+RadiusNoseRotorOutlet^2)/2)^0.5; %m
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
    ((MassFlowParameterRotor./(2*pi.*RadiusRotor.*WidthRotorInlet))*ones(size(MachNoRotorInletGuess))).*(1+((GammaTurbine-1)./2)*MachNoRotorInletGuess.^2).^((1./(GammaTurbine-1))*ones(size(MachNoRotorInletGuess)));
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
VelocityMeridionalRotorInlet = MassFlowParameterRotor./(2*pi.*RadiusRotor.*WidthRotorInlet).*(1+(GammaTurbine-1)./2.*MachNoRotorInlet.^2).^(1./(GammaTurbine-1)); %m/sec
AngleRelativeRotorInlet = atan(VelocityTangentialRelativeRotorInlet./VelocityMeridionalRotorInlet); %rad
FlowCoefficient = VelocityMeridionalRotorInlet./SpeedRotorTip;
%
%section 2 velocity triangle
TemperatureRatioTurbine = 1-EnthalpySpecificTurbine./CpTurbine./TemperatureTurbineInlet;
VelocityMeridionalRotorOutlet = (MassFlowParameterRotor./AreaRotorOutlet./TemperatureRatioTurbine.*PressureRatioTurbine); %m/sec
% VelocityTangentialRotorOutlet = SpeedRotorTip.*(RadiusRMSRotorOutlet/RadiusRotor)-VelocityMeridionalRotorOutlet.*tan(AngleBladeRotorOutlet); %m/sec
VelocityRelativeRotorOutlet = VelocityMeridionalRotorOutlet./cos(AngleBladeRotorOutlet); %m/sec
%find rotor outlet relative Mach No
TemperatureTurbineOutlet = TemperatureTurbineInlet./TemperatureRatioTurbine; %K
MachNoRelativeRotorOutlet = VelocityRelativeRotorOutlet./(GammaTurbineOutlet.*RTurbine.*TemperatureTurbineOutlet).^0.5;
%
%
%% fit peak efficiency incidence angle vs. tip speed
%find data from different training maps
% MapBreakpoints = 1;
% for k = 2:length(AreaRadiusRatioMultiplier)
%     if abs(AreaRadiusRatioMultiplier(k)-AreaRadiusRatioMultiplier(k-1)) > 0.04 %4% change in A/R ratio
%         MapBreakpoints = [MapBreakpoints,k];
%     end
% end
% MapBreakpoints = [MapBreakpoints,length(AreaRadiusRatioMultiplier)+1];
% %find speedlines and fit peak efficiency for each data set
% CoeffsAngleIncidenceVsSpeedEfficiencyPeak = zeros(2,length(MapBreakpoints)-1);
% CoeffsFlowCoefficientVsSpeedEfficiencyPeak = zeros(2,length(MapBreakpoints)-1);
% MachNoChoke = zeros(1,length(MapBreakpoints)-1);
% SpeedRotorTipEfficiencyMax = zeros(1,length(MapBreakpoints)-1);
% for k1 = 1:length(MapBreakpoints)-1
    %find speed breakpoints
    SpeedRotorTurbineMap = SpeedRotorTip;%(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
    EfficiencyTurbineMap = EfficiencyTurbine;%(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
    AngleRelativeRotorInletMap = AngleRelativeRotorInlet;%(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
    FlowCoefficientMap = FlowCoefficient;%(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
    MachNoRelativeRotorOutletMap = MachNoRelativeRotorOutlet;%(MapBreakpoints(k1):MapBreakpoints(k1+1)-1);
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
    
%correct optimum rotor inlet incidence for A/R    
AngleIncidenceEfficiencyPeak = atan(AreaRadiusRatioMultiplier(MaxEfficiencyIndex).*(1./FlowCoefficientMap(MaxEfficiencyIndex)+...
    tan(AngleRelativeRotorInletMap(MaxEfficiencyIndex)))-1./FlowCoefficientMap(MaxEfficiencyIndex));
AngleIncidenceEfficiencyPeak = (AngleIncidenceEfficiencyPeak+AngleRelativeRotorInletMap(MaxEfficiencyIndex))./2;
%fit peak efficiency incidence angle vs. tip speed
CoeffsAngleIncidenceVsSpeedEfficiencyPeak(:,1) = [SpeedRotorTurbineMap(MaxEfficiencyIndex),ones(length(MaxEfficiencyIndex),1)]\AngleIncidenceEfficiencyPeak;
CoeffsAngleIncidenceVsSpeedEfficiencyPeak(:,2) = [SpeedRotorTurbineMap(MaxEfficiencyIndex),ones(length(MaxEfficiencyIndex),1)]\FlowCoefficientMap(MaxEfficiencyIndex);
% AngleIncidencePeakModel = [SpeedRotorTurbineMap(MaxEfficiencyIndex),ones(length(MaxEfficiencyIndex),1)]*CoeffsAngleIncidenceVsSpeedEfficiencyPeak;
% figure(1)
% plot(SpeedRotorTurbineMap(MaxEfficiencyIndex),AngleIncidenceEfficiencyPeak,'o',SpeedRotorTurbineMap(MaxEfficiencyIndex),AngleIncidencePeakModel(:,1))
% figure(2)
% plot(SpeedRotorTurbineMap(MaxEfficiencyIndex),FlowCoefficientMap(MaxEfficiencyIndex),'o',SpeedRotorTurbineMap(MaxEfficiencyIndex),AngleIncidencePeakModel(:,2))
%
%
%% fit max efficiency rotor speed
    MachNoChoke = max(MachNoRelativeRotorOutletMap);
    MaxEfficiency = max(EfficiencyTurbineMap(MaxEfficiencyIndex));
    for k = 1:length(MaxEfficiencyIndex)
        if EfficiencyTurbineMap(MaxEfficiencyIndex(k)) == MaxEfficiency
            SpeedRotorTipEfficiencyMax = SpeedRotorTurbineMap(MaxEfficiencyIndex(k));
        end
    end
%end
%
%
%
%
end