function TurbineModel = TurbineMap7(NumberTrainingMaps,NumberValidationMaps)
%
%
%
%
%% initialize class
TurbineModel = TurbineModel7;
%
%
%% load training data
for k = 1:NumberTrainingMaps
    fprintf('Select Training Turbine\n')
    [FileNameTurbine,PathNameTurbine] = uigetfile('*.*','Select Turbine'); %prompt to get turbine raw data
    TurbineModel = TurbineModel.TurbineMapDataPrep([PathNameTurbine,FileNameTurbine],'Training');
end
%
%
%% load validation data
if NumberValidationMaps > 0
    for k = 1:NumberValidationMaps
        fprintf('Select Validation Turbine\n')
        [FileNameTurbine,PathNameTurbine] = uigetfile('*.*','Select Turbine'); %prompt to get turbine raw data
        TurbineModel = TurbineModel.TurbineMapDataPrep([PathNameTurbine,FileNameTurbine],'Validation');
    end
end
%
%
%% build turbine model
TurbineModel = TurbineModel.BuildTurbineModel;
TurbineModel = TurbineModel.SolveTurbineModel; %extend speedlines
%
%
%% plot models
TurbineModel.PlotModel
%
%
end