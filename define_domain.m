%% Define domain
% computes important variables that describe the domain: variable radius and its derivative,
% collocation points, normal vectors, source points, area and center of mass

function [r, dr_dt, Gamma_col, Normal_col, p_source, A_Gamma, CoM] = define_domain(a,b,t_col,p_beta)

% variable radius
[r, dr_dt] = define_r(a,b);
       
% collocation points on domain boundary
Gamma_col = [r(t_col) .* cos(t_col); r(t_col) .* sin(t_col)]; 

% normal vectors in collocation points
Normal_temp = [ dr_dt(t_col) .* sin(t_col) + r(t_col) .* cos(t_col); -dr_dt(t_col) .* cos(t_col) + r(t_col) .* sin(t_col)]; 
Normal_norm = sqrt(Normal_temp(1,:).^2+ Normal_temp(2,:).^2);
Normal_col = [Normal_temp(1,:)./Normal_norm ; Normal_temp(2,:)./Normal_norm ];

% source points
p_source = Gamma_col + p_beta * Normal_col;

% compute are and center of mass
dt_col = t_col(2)-t_col(1);
A_Gamma = dt_col/2*(sum(r(t_col).^2)); 
CoM = dt_col/(3*A_Gamma)*[sum(cos(t_col).*r(t_col).^3); sum(sin(t_col).*r(t_col).^3)];
    
% % plot domain and point sources
% figure
% plot(Gamma_col(1,:),Gamma_col(2,:))
% hold on
% plot(p_source(1,:),p_source(2,:),'r')
% hold off
% grid on

end