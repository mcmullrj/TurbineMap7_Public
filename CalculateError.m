function ErrorTraining = CalculateError(TrainingData,TurbineDataModelAll,ConvergenceIndexTraining)
%turbine model error calculations
%created by Bob McMullen
%
%
%calculate training data error
ErrorMassFlowTraining = TurbineDataModelAll(:,2)-TrainingData(:,3);
ErrorEfficiencyTraining = (TurbineDataModelAll(:,4)-TrainingData(:,5));
ErrorAllTraining = [ErrorMassFlowTraining./TrainingData(:,3).*100,ErrorEfficiencyTraining./TrainingData(:,5).*100];
ErrorStdDevTraining = std(ErrorAllTraining);
ErrorMeanAllTraining = (ones(1,length(ConvergenceIndexTraining))*ErrorAllTraining(ConvergenceIndexTraining,:))./length(ConvergenceIndexTraining);
ErrorRMSAllTraining = ((ones(1,length(ConvergenceIndexTraining))*ErrorAllTraining(ConvergenceIndexTraining,:).^2)./length(ConvergenceIndexTraining)).^0.5;
UnchokeIndex = [];
ChokeIndex = [];   
for k = 1:length(ConvergenceIndexTraining)
    if TurbineDataModelAll(ConvergenceIndexTraining(k),6) == 0
        UnchokeIndex = [UnchokeIndex;ConvergenceIndexTraining(k)];
    else
        ChokeIndex = [ChokeIndex;ConvergenceIndexTraining(k)];
    end
end
MaxMassFlow = max(TrainingData(:,3));
LowFlowIndex = [];
UnchokeIndex = [];
for k = 1:length(ConvergenceIndexTraining)
    if TrainingData(ConvergenceIndexTraining(k),3) < 0.9*MaxMassFlow
        LowFlowIndex = [LowFlowIndex;k];
    elseif TurbineDataModelAll(ConvergenceIndexTraining(k),6) == 0
        UnchokeIndex = [UnchokeIndex;k];
    end
end
ErrorRMSUnchokeTraining = ((ones(1,length(UnchokeIndex))*ErrorAllTraining(UnchokeIndex,:).^2)./length(UnchokeIndex)).^0.5;
ErrorRMSChokeTraining = ((ones(1,length(ChokeIndex))*ErrorAllTraining(ChokeIndex,:).^2)./length(ChokeIndex)).^0.5;
ErrorMeanUnchokeTraining = (ones(1,length(UnchokeIndex))*ErrorAllTraining(UnchokeIndex,:))./length(UnchokeIndex);
ErrorStdDevUnchokeTraining = std(ErrorAllTraining(UnchokeIndex,:));
ErrorRMSLowFlowTraining = ((ones(1,length(LowFlowIndex))*ErrorAllTraining(LowFlowIndex,:).^2)./length(LowFlowIndex)).^0.5;
ErrorMeanLowFlowTraining = (ones(1,length(LowFlowIndex))*ErrorAllTraining(LowFlowIndex,:))./length(LowFlowIndex);
ErrorStdDevLowFlowTraining = std(ErrorAllTraining(LowFlowIndex,:));        
ErrorTraining = [ErrorMeanAllTraining;ErrorStdDevTraining;ErrorRMSAllTraining;...
    ErrorMeanUnchokeTraining;ErrorStdDevUnchokeTraining;ErrorRMSUnchokeTraining;...
    ErrorMeanLowFlowTraining;ErrorStdDevLowFlowTraining;ErrorRMSLowFlowTraining]';
%
%
%
%
end