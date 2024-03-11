function Error2MassFlowCorrected = FitTurbineGeometry7(DesignVariables)
%geometry optimization function
%created by Bob McMullen
%
%
%% design variables
AreaRadiusRatioVoluteThroat = DesignVariables(1); %m
RadiusNoseRotorOutlet = DesignVariables(2); %m
AngleBladeRotorOutlet = DesignVariables(3); %rad
%
%
%% turbine map and geometry inputs
global GeometryTurbine TurbineDataVPSortGoodData
%extract geometry
RadiusRotor = GeometryTurbine(:,1);
RadiusTipRotorOutlet = GeometryTurbine(:,2);
AreaRadiusRatioScaleTraining = GeometryTurbine(:,4);
%extract turbine performance data
SpeedRotorTip = TurbineDataVPSortGoodData(:,8); %m/sec
MassFlowCorrected = TurbineDataVPSortGoodData(:,3)./1000; %kg/sec*sqrtK/Paa
PressureRatioTurbine = TurbineDataVPSortGoodData(:,4);
EfficiencyTurbine = TurbineDataVPSortGoodData(:,5);
TemperatureTurbineInlet = TurbineDataVPSortGoodData(:,6); %K
RTurbine = TurbineDataVPSortGoodData(:,10); %J/kg/K
CpTurbine = TurbineDataVPSortGoodData(:,9); %J/kg/K
GammaTurbine = TurbineDataVPSortGoodData(:,11);
%
%
%% update turbine geometry
RadiusRMSRotorOutlet = ((RadiusTipRotorOutlet.^2+RadiusNoseRotorOutlet^2)./2).^0.5; %m
AreaRotorOutlet = pi.*(RadiusTipRotorOutlet.^2-RadiusNoseRotorOutlet^2); %m^2
%
%
%% Euler turbomachine equation
MassFlowCorrectedModel = ((CpTurbine.*TemperatureTurbineInlet.*EfficiencyTurbine./SpeedRotorTip.*(1-PressureRatioTurbine.^((1-GammaTurbine)./GammaTurbine))+(RadiusRMSRotorOutlet./RadiusRotor).^2.*SpeedRotorTip)./...
    (1./RadiusRotor./(AreaRadiusRatioScaleTraining.*AreaRadiusRatioVoluteThroat)+RadiusRMSRotorOutlet./RadiusRotor.*tan(AngleBladeRotorOutlet)./AreaRotorOutlet.*(1-EfficiencyTurbine.*(1-PressureRatioTurbine.^((1-GammaTurbine)./GammaTurbine))).*PressureRatioTurbine))./...
    RTurbine./TemperatureTurbineInlet.^0.5;
%
%
%% calculate corrected flow error
MaxMassFlowCorrected = max(MassFlowCorrected);
ErrorMassFlowCorrected = (MassFlowCorrectedModel-MassFlowCorrected)./MaxMassFlowCorrected; %normalize error before squaring to avoid error^2 meeting convergence tolerance
% figure(2)
% plot(PressureRatioTurbine,ErrorMassFlowCorrected,'o')
Error2MassFlowCorrected = ((ones(1,length(ErrorMassFlowCorrected))*ErrorMassFlowCorrected.^2)/length(ErrorMassFlowCorrected));
%
%
%
%
end