%% Define r
% defines the variable radius and its derivative for given Fourier
% coefficients

function [r, dr_dt] = define_r(a,b)

r =  @(t) (cos(kron((0:length(a)-1),t'))*a' + sin(kron((0:length(b)-1),t'))*b')'; 
dr_dt = @(t) (-sin(kron((0:length(a)-1),t'))*((0:length(a)-1).*a)' + cos(kron((0:length(b)-1),t'))*((0:length(b)-1).*b)')';

end