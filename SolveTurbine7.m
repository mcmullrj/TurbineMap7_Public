function TurbinePrediction = SolveTurbine7(GeometryTurbine,CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsLossModel,MachNoChoke,SpeedRotorTipEfficiencyMax,TurbineDataInput,ConvergenceEfficiencyError)
%turbine model solver
%created by Bob McMullen
%
%
%% design variables
AreaRadiusRatioMultiplier = GeometryTurbine(:,4);
AngleBladeRotorOutlet = GeometryTurbine(:,6); %rad
RadiusRotor = GeometryTurbine(:,1); %m
RadiusTipRotorOutlet = GeometryTurbine(:,2); %m
RadiusNoseRotorOutlet = GeometryTurbine(:,5); %m
WidthRotorInlet = GeometryTurbine(:,7); %m
AreaRadiusRatioVoluteBaseModel = GeometryTurbine(:,3); %m^2
AreaRadiusRatioVoluteThroat = AreaRadiusRatioMultiplier.*AreaRadiusRatioVoluteBaseModel; %m^2
%
%
%% turbine map and geometry inputs
%extract turbine performance data
SpeedRotorTip = TurbineDataInput(:,1); %m/sec
PressureRatioTurbine = TurbineDataInput(:,2);
TemperatureTurbineInlet = TurbineDataInput(:,3); %K
RTurbine = TurbineDataInput(:,5); %J/kg/K
CpTurbine = TurbineDataInput(:,4); %J/kg/K
GammaTurbine = TurbineDataInput(:,6);
%
%
%% update turbine geometry
RadiusRMSRotorOutlet = ((RadiusTipRotorOutlet.^2+RadiusNoseRotorOutlet.^2)./2).^0.5; %m
RadiusRotorPassage = RadiusTipRotorOutlet-RadiusNoseRotorOutlet; %m
AreaRotorOutlet = pi.*(RadiusTipRotorOutlet.^2-RadiusNoseRotorOutlet.^2); %m^2
% AngleIncidenceEfficiencyPeak = [SpeedRotorTip.^2,SpeedRotorTip,ones(size(SpeedRotorTip))]*CoeffsAngleIncidenceVsSpeedEfficiencyPeak;
% AngleIncidenceEfficiencyPeak = ones(size(SpeedRotorTip)).*CoeffsAngleIncidenceVsSpeedEfficiencyPeak;
%
%
%% solve unchoked flow model - Newton's method for turbine efficiency error
EfficiencyTurbineInitialGuessHigh = 0.9;
EfficiencyTurbineInitialGuessLow = 0.2; %add a second search from low efficiency?
ChokeFlag = 0;
dEfficiency = 1e-7;
IterationsMax = 100;
%start search from high efficiency
EfficiencyTurbineGuessAll = EfficiencyTurbineInitialGuessHigh+[-2*dEfficiency,-dEfficiency,0,dEfficiency,2*dEfficiency];
ErrorEfficiencyAll = zeros(length(SpeedRotorTip),length(EfficiencyTurbineGuessAll));
%evaluate initial efficiency error and first and second derivative
for k1 = 1:length(EfficiencyTurbineGuessAll)
    [ErrorEfficiencyAll(:,k1),~] = TurbineLossModel(SpeedRotorTip,PressureRatioTurbine,EfficiencyTurbineGuessAll(k1),TemperatureTurbineInlet,CpTurbine,GammaTurbine,RTurbine,...
        RadiusRotor,RadiusRMSRotorOutlet,RadiusRotorPassage,AreaRadiusRatioVoluteThroat,AreaRadiusRatioVoluteBaseModel,AreaRotorOutlet,AngleBladeRotorOutlet,WidthRotorInlet,CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsLossModel,MachNoChoke,ChokeFlag,SpeedRotorTipEfficiencyMax);
end
%determine efficiency step
dErrorEfficiency = (ErrorEfficiencyAll(:,4)-ErrorEfficiencyAll(:,2))./(2*dEfficiency);
d2ErrorEfficiency = (ErrorEfficiencyAll(:,5)-2.*ErrorEfficiencyAll(:,3)+ErrorEfficiencyAll(:,1))./(4*dEfficiency^2);
EfficiencyTurbineGuessAllUpdate = ones(size(SpeedRotorTip))*EfficiencyTurbineGuessAll(3)-dErrorEfficiency./d2ErrorEfficiency;
%solve each data point
EfficiencyTurbineModel = zeros(size(SpeedRotorTip));
for k2 = 1:length(SpeedRotorTip)
    EfficiencyTurbineGuessUpdate = EfficiencyTurbineGuessAllUpdate(k2);
    ErrorEfficiency = ErrorEfficiencyAll(k2,:);
    EfficiencyTurbineGuess = EfficiencyTurbineGuessAll;
    IterationCount = 1;
    while (IterationCount < IterationsMax) && (ErrorEfficiency(3) > 0.1*ConvergenceEfficiencyError) && (abs(EfficiencyTurbineGuessUpdate-EfficiencyTurbineGuess(3)) > ConvergenceEfficiencyError)
        EfficiencyTurbineGuess = EfficiencyTurbineGuessUpdate+[-2*dEfficiency,-dEfficiency,0,dEfficiency,2*dEfficiency];
        for k1 = 1:length(EfficiencyTurbineGuess)
            [ErrorEfficiency(k1),~] = TurbineLossModel(SpeedRotorTip(k2),PressureRatioTurbine(k2),EfficiencyTurbineGuess(k1),TemperatureTurbineInlet(k2),CpTurbine(k2),GammaTurbine(k2),RTurbine(k2),...
                RadiusRotor(k2),RadiusRMSRotorOutlet(k2),RadiusRotorPassage(k2),AreaRadiusRatioVoluteThroat(k2),AreaRadiusRatioVoluteBaseModel(k2),AreaRotorOutlet(k2),AngleBladeRotorOutlet(k2),WidthRotorInlet(k2),CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsLossModel,MachNoChoke,ChokeFlag,SpeedRotorTipEfficiencyMax);
        end
        dErrorEfficiency = (ErrorEfficiency(4)-ErrorEfficiency(2))/(2*dEfficiency);
        d2ErrorEfficiency = (ErrorEfficiency(5)-2*ErrorEfficiency(3)+ErrorEfficiency(1))/(4*dEfficiency^2);
        EfficiencyTurbineGuessUpdate = EfficiencyTurbineGuess(3)-dErrorEfficiency/d2ErrorEfficiency;
        IterationCount = IterationCount+1;
    end
    EfficiencyTurbineModel(k2) = EfficiencyTurbineGuessUpdate;
end
%get final solution
[ErrorEfficiencyFinal,MassFlowCorrectedModel,ChokedFlowIndex] = TurbineLossModel(SpeedRotorTip,PressureRatioTurbine,EfficiencyTurbineModel,TemperatureTurbineInlet,CpTurbine,GammaTurbine,RTurbine,...
    RadiusRotor,RadiusRMSRotorOutlet,RadiusRotorPassage,AreaRadiusRatioVoluteThroat,AreaRadiusRatioVoluteBaseModel,AreaRotorOutlet,AngleBladeRotorOutlet,WidthRotorInlet,CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsLossModel,MachNoChoke,ChokeFlag,SpeedRotorTipEfficiencyMax);
TurbineDataOutput = [SpeedRotorTip,MassFlowCorrectedModel.*1000,PressureRatioTurbine,EfficiencyTurbineModel,ErrorEfficiencyFinal,zeros(size(SpeedRotorTip))];
TurbineDataOutput(ChokedFlowIndex,6) = 1;
%
%
%% solve choked flow data points efficiency
% ChokedFlowIndex = [];
if ~isempty(ChokedFlowIndex)
    ChokeFlag = 1;
    ErrorEfficiencyAll = zeros(length(ChokedFlowIndex),length(EfficiencyTurbineGuessAll));
    %evaluate initial efficiency error and first and second derivative
    for k1 = 1:length(EfficiencyTurbineGuessAll)
        [ErrorEfficiencyAll(:,k1),~] = TurbineLossModel(SpeedRotorTip(ChokedFlowIndex),PressureRatioTurbine(ChokedFlowIndex),EfficiencyTurbineGuessAll(k1),TemperatureTurbineInlet(ChokedFlowIndex),CpTurbine(ChokedFlowIndex),GammaTurbine(ChokedFlowIndex),RTurbine(ChokedFlowIndex),...
            RadiusRotor(ChokedFlowIndex),RadiusRMSRotorOutlet(ChokedFlowIndex),RadiusRotorPassage(ChokedFlowIndex),AreaRadiusRatioVoluteThroat(ChokedFlowIndex),AreaRadiusRatioVoluteBaseModel(ChokedFlowIndex),AreaRotorOutlet(ChokedFlowIndex),AngleBladeRotorOutlet(ChokedFlowIndex),WidthRotorInlet(ChokedFlowIndex),...
            CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsLossModel,MachNoChoke,ChokeFlag,SpeedRotorTipEfficiencyMax);
    end
    %determine efficiency step
    dErrorEfficiency = (ErrorEfficiencyAll(:,4)-ErrorEfficiencyAll(:,2))./(2*dEfficiency);
    d2ErrorEfficiency = (ErrorEfficiencyAll(:,5)-2.*ErrorEfficiencyAll(:,3)+ErrorEfficiencyAll(:,1))./(4*dEfficiency^2);
    EfficiencyTurbineGuessAllUpdate = ones(length(ChokedFlowIndex),1)*EfficiencyTurbineGuessAll(3)-dErrorEfficiency./d2ErrorEfficiency;
    %solve each data point
    EfficiencyTurbineModel = zeros(length(ChokedFlowIndex),1);
    for k2 = 1:length(ChokedFlowIndex)
        EfficiencyTurbineGuessUpdate = EfficiencyTurbineGuessAllUpdate(k2);
        ErrorEfficiency = ErrorEfficiencyAll(k2,:);
        EfficiencyTurbineGuess = EfficiencyTurbineGuessAll;
        IterationCount = 1;
        while (IterationCount < IterationsMax) && (ErrorEfficiency(3) > 0.1*ConvergenceEfficiencyError) && (abs(EfficiencyTurbineGuessUpdate-EfficiencyTurbineGuess(3)) > ConvergenceEfficiencyError)
            EfficiencyTurbineGuess = EfficiencyTurbineGuessUpdate+[-2*dEfficiency,-dEfficiency,0,dEfficiency,2*dEfficiency];
            for k1 = 1:length(EfficiencyTurbineGuess)
                [ErrorEfficiency(k1),~] = TurbineLossModel(SpeedRotorTip(ChokedFlowIndex(k2)),PressureRatioTurbine(ChokedFlowIndex(k2)),EfficiencyTurbineGuess(k1),TemperatureTurbineInlet(ChokedFlowIndex(k2)),CpTurbine(ChokedFlowIndex(k2)),GammaTurbine(ChokedFlowIndex(k2)),RTurbine(ChokedFlowIndex(k2)),...
                    RadiusRotor(ChokedFlowIndex(k2)),RadiusRMSRotorOutlet(ChokedFlowIndex(k2)),RadiusRotorPassage(ChokedFlowIndex(k2)),AreaRadiusRatioVoluteThroat(ChokedFlowIndex(k2)),AreaRadiusRatioVoluteBaseModel(ChokedFlowIndex(k2)),AreaRotorOutlet(ChokedFlowIndex(k2)),AngleBladeRotorOutlet(ChokedFlowIndex(k2)),WidthRotorInlet(ChokedFlowIndex(k2)),...
                    CoeffsAngleIncidenceVsSpeedEfficiencyPeak,CoeffsLossModel,MachNoChoke,ChokeFlag,SpeedRotorTipEfficiencyMax);
            end
            dErrorEfficiency = (ErrorEfficiency(4)-ErrorEfficiency(2))/(2*dEfficiency);
            d2ErrorEfficiency = (ErrorEfficiency(5)-2*ErrorEfficiency(3)+ErrorEfficiency(1))/(4*dEfficiency^2);
            EfficiencyTurbineGuessUpdate = EfficiencyTurbineGuess(3)-dErrorEfficiency/d2ErrorEfficiency;
            IterationCount = IterationCount+1;
        end
        EfficiencyTurbineModel(k2) = EfficiencyTurbineGuessUpdate;
    end
    TurbineDataOutput(ChokedFlowIndex,4) = EfficiencyTurbineModel;
end
TurbinePrediction = TurbineDataOutput;
%
%
end %SolveTurbine7
%
%
%
%
function [ErrorEfficiency,MassFlowCorrectedModel,ChokedFlowIndex] = TurbineLossModel(SpeedRotorTip,PressureRatioTurbine,EfficiencyTurbineGuess,TemperatureTurbineInlet,CpTurbine,GammaTurbine,RTurbine,...
    RadiusRotor,RadiusRMSRotorOutlet,RadiusRotorPassage,AreaRadiusRatioVoluteThroat,AreaRadiusRatioVoluteBaseModel,AreaRotorOutlet,AngleBladeRotorOutlet,WidthRotorInlet,CoeffsAngleIncidenceVsSpeedEfficiencyPeak,...
    CoeffsEnthalpyLossModel,MachNoChoke,ChokeFlag,SpeedRotorTipEfficiencyMax)
%
%
EnthalpySpecificTurbineIdeal = CpTurbine.*TemperatureTurbineInlet.*(1-PressureRatioTurbine.^((1-GammaTurbine)./GammaTurbine)); %J/kg
%% Euler turbomachine equation and velocity triangles
MassFlowCorrectedModel = ((CpTurbine.*TemperatureTurbineInlet.*EfficiencyTurbineGuess./SpeedRotorTip.*(1-PressureRatioTurbine.^((1-GammaTurbine)./GammaTurbine))+(RadiusRMSRotorOutlet./RadiusRotor).^2.*SpeedRotorTip)./...
    (1./RadiusRotor./AreaRadiusRatioVoluteThroat+RadiusRMSRotorOutlet./RadiusRotor.*tan(AngleBladeRotorOutlet)./AreaRotorOutlet.*(1-EfficiencyTurbineGuess.*(1-PressureRatioTurbine.^((1-GammaTurbine)./GammaTurbine))).*PressureRatioTurbine))./...
    RTurbine./TemperatureTurbineInlet.^0.5;
% section 1 velocity triangle
MassFlowParameterRotor = MassFlowCorrectedModel.*RTurbine.*TemperatureTurbineInlet.^0.5;
VelocityTangentialRotorInlet = MassFlowParameterRotor./AreaRadiusRatioVoluteThroat./RadiusRotor; %m/sec
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
VelocityRotorInlet = (VelocityTangentialRotorInlet.^2+VelocityMeridionalRotorInlet.^2).^0.5;
VelocityRelativeRotorInlet = (VelocityTangentialRelativeRotorInlet.^2+VelocityMeridionalRotorInlet.^2).^0.5; %m/sec
AngleRelativeRotorInlet = atan(VelocityTangentialRelativeRotorInlet./VelocityMeridionalRotorInlet); %rad
%correct peak efficiency points for incidence
AngleIncidenceEfficiencyPeakBase = [SpeedRotorTip,ones(size(SpeedRotorTip))]*CoeffsAngleIncidenceVsSpeedEfficiencyPeak; %column 1 = optimum angle (rad), column 2 = optimum flow coefficient
AngleIncidenceEfficiencyPeak = atan(AreaRadiusRatioVoluteBaseModel./AreaRadiusRatioVoluteThroat.*(1./AngleIncidenceEfficiencyPeakBase(:,2)+tan(AngleIncidenceEfficiencyPeakBase(:,1)))-1./AngleIncidenceEfficiencyPeakBase(:,2));
AngleIncidenceEfficiencyPeak = (AngleIncidenceEfficiencyPeak+AngleIncidenceEfficiencyPeakBase(:,1))./2;
% EfficiencyRatioPeak = (1+AngleIncidenceEfficiencyPeakBase(:,2).*tan(AngleIncidenceEfficiencyPeak))./(1+AngleIncidenceEfficiencyPeakBase(:,2).*tan(AngleIncidenceEfficiencyPeakBase(:,1)));
%section 2 velocity triangle
EnthalpySpecificTurbine = EfficiencyTurbineGuess.*EnthalpySpecificTurbineIdeal; %J/kg
TemperatureRatioTurbine = 1-EnthalpySpecificTurbine./CpTurbine./TemperatureTurbineInlet;
VelocityMeridionalRotorOutlet = (MassFlowParameterRotor./AreaRotorOutlet./TemperatureRatioTurbine.*PressureRatioTurbine); %m/sec
VelocityRelativeRotorOutlet = VelocityMeridionalRotorOutlet./cos(AngleBladeRotorOutlet); %m/sec
%
%
%% modify rotor choked flow conditions
TemperatureTurbineOutlet = TemperatureTurbineInlet./TemperatureRatioTurbine; %K
%update - calc turbine out Cp, gamma
GammaTurbineOutlet = GammaTurbine;
TemperatureRelativeRotorOutlet = VelocityRelativeRotorOutlet./(RTurbine.*GammaTurbineOutlet.*MachNoChoke^2);
VelocityRelativeRotorOutletChoke = MachNoChoke.*(GammaTurbineOutlet.*RTurbine.*TemperatureTurbineOutlet).^0.5; %m/sec
VelocityMeridionalRotorOutletMax = VelocityRelativeRotorOutletChoke.*cos(AngleBladeRotorOutlet);
PressureRatioRotorOutlet = MassFlowParameterRotor./AreaRotorOutlet./TemperatureRatioTurbine.*PressureRatioTurbine./VelocityMeridionalRotorOutletMax;
EfficiencyLossChoke = zeros(size(SpeedRotorTip));
ChokedFlowIndex = [];
for k = 1:length(SpeedRotorTip)
    if VelocityMeridionalRotorOutlet(k) > VelocityMeridionalRotorOutletMax(k)
        ChokedFlowIndex = [ChokedFlowIndex;k];
        if ChokeFlag == 1
            EfficiencyLossChoke(k) = (1-PressureRatioRotorOutlet(k)^((1-GammaTurbine(k))/GammaTurbine(k)))/(1-PressureRatioTurbine(k)^((1-GammaTurbine(k))/GammaTurbine(k)));
%             EfficiencyLossChoke(k) = CpTurbine(k).*(TemperatureRelativeRotorOutlet(k)-TemperatureTurbineOutlet(k))./EnthalpySpecificTurbineIdeal(k);
        end
    end
end
%
%
%% enthalpy losses
EnthalpyLossModelIncidence = (sin(AngleRelativeRotorInlet-AngleIncidenceEfficiencyPeak).^2.*VelocityRelativeRotorInlet.^2./2)./EnthalpySpecificTurbineIdeal; %J/kg
EnthalpyLossModelRotorFriction = (cos(AngleRelativeRotorInlet-AngleIncidenceEfficiencyPeak).^2.*VelocityRelativeRotorInlet.^2./2+VelocityRelativeRotorOutlet.^2./2)./EnthalpySpecificTurbineIdeal; %J/kg
% EnthalpyLossModelOutletSwirl = (VelocityTangentialRotorOutlet.^2./2)./EnthalpySpecificTurbineIdeal;
% EnthalpyLossModelBearings = SpeedRotorTip.^2./(1e5.*MassFlowCorrectedModel.*(PressureRatioTurbine.*101300)./TemperatureTurbineInlet.^0.5)./EnthalpySpecificTurbineIdeal;
EnthalpyLossModelStator = VelocityRotorInlet.^2./2./EnthalpySpecificTurbineIdeal;
EnthalpyLossPassageVortex = 10.*0.5.*RadiusRotorPassage.*...
    ((SpeedRotorTipEfficiencyMax-SpeedRotorTip)./SpeedRotorTipEfficiencyMax)./(pi.*RadiusRotor.^2).*MassFlowCorrectedModel.*RTurbine.*TemperatureTurbineInlet.^0.5.*SpeedRotorTip.*...
    (1./AreaRadiusRatioVoluteThroat-1./AreaRadiusRatioVoluteBaseModel)./EnthalpySpecificTurbineIdeal; %default calculation is for increasing A/R
EnthalpyLossPassageVortexDecreasing = 10.*0.5.*RadiusRotorPassage.*...
    ((SpeedRotorTipEfficiencyMax-SpeedRotorTip)./SpeedRotorTip)./(pi.*RadiusRotor.^2).*MassFlowCorrectedModel.*RTurbine.*TemperatureTurbineInlet.^0.5.*SpeedRotorTip.*...
    (1./AreaRadiusRatioVoluteThroat-1./AreaRadiusRatioVoluteBaseModel)./EnthalpySpecificTurbineIdeal;
for k = 1:length(SpeedRotorTip)
    if AreaRadiusRatioVoluteThroat(k) < AreaRadiusRatioVoluteBaseModel(k)
        EnthalpyLossPassageVortex(k) = EnthalpyLossPassageVortexDecreasing(k);
    end
end
%
%
%% efficiency model
% EnthalpyLossModel = [EnthalpyLossModelBearings,EnthalpyLossModelIncidence,EnthalpyLossModelRotorFriction,EnthalpyLossModelOutletSwirl,ones(size(EnthalpyLossModelBearings))]*CoeffsEnthalpyLossModel; %J/kg
EnthalpyLossModel = [EnthalpyLossModelIncidence,EnthalpyLossModelRotorFriction,EnthalpyLossModelStator]*CoeffsEnthalpyLossModel; %J/kg
% EnthalpyLossModel = [EnthalpyLossModelBearings,EnthalpyLossModelIncidence,EnthalpyLossModelRotorFriction,EnthalpyLossModelStator]*CoeffsEnthalpyLossModel; %J/kg
EfficiencyModel = (1-EnthalpyLossModel-EnthalpyLossPassageVortex-EfficiencyLossChoke);
ErrorEfficiency = (EfficiencyTurbineGuess-EfficiencyModel).^2;
%
%
end