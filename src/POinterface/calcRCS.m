function [RCS] = calcRCS(Const, Solver_setup, theta_grid, phi_grid, xVectors)

%calcRCS
%   Usage:
%       [RCS] = calcRCS(Const, Solver_setup, theta_grid, phi_grid, xVectors)
%   Input Arguments:
%       Const
%           A global struct, containing general data
%       Solver_setup
%           Solver specific struct, e.g. frequency range, basis function details, geometry details
%       theta_grid
%           List of theta incident angles
%       phi_grid
%           List of phi incident angles
%       xVectors
%           The reference solution-vector data (e.g. MoM solution of FEKO or SUN-EM)
%
%   Output Arguments:
%       RCS
%           Structs containing monostatic RCS solution and timing data
%
%   Description:
%       Runs the PO solution based on the RGB basis functions parsed
%       from FEKO for each incident angle and calculates the RCS
%
%   =======================
%   Written by Cullen Stewart-Burger on 27 July 2020
%   Stellenbosch University
%   Email: 20751028@sun.ac.za

narginchk(5,5);
r = 100000;
Solver_setup.r_bounding = getBoundingRadius(Solver_setup);
% Calculate now the E-field value here internal
index = 0;
RCS = zeros(length(theta_grid)*length(phi_grid), 1);
Raytracer.setGeom(Solver_setup);
Solver_setup.Visibility_matrix = selfShadow(Solver_setup);

for theta_degrees = theta_grid
    for phi_degrees = phi_grid
        index = index + 1;
        
        Solver_setup.theta = theta_degrees;
        Solver_setup.phi = phi_degrees;
        
        % --------------------------------------------------------------------------------------------------
        % Run the EM solver
        % --------------------------------------------------------------------------------------------------
        %         tic;
        %         [Vpp Vpn Vnp Vnn] = selfShadow(Solver_setup);
        %         toc;
        
        [Solution] = runEMsolvers(Const, Solver_setup, 0, 0, xVectors);
        EfieldAtPointSpherical =  calculateEfieldAtPointRWG(Const, r, theta_degrees, phi_degrees, ...
            Solver_setup, Solution.PO.Isol);
        
%                 HfieldAtPointSpherical = zeros(10/0.01+1, 3);
%                 r_obs = 0:0.001:1;
%                 for ind = 1:1001
%                     HfieldAtPointSpherical(ind, :) =  calculateHfieldAtPointRWGCart(Const, [0, 0, r_obs(ind)], ...
%                         Solver_setup, Solution.PO.Isol);
% %                     HfieldAtPointSpherical(ind, :) =  calculateHfieldAtPointRWG(Const, r_obs(ind), 0, 0, ...
% %                         Solver_setup, Solution.PO.Isol);
%                 end
%         
%                 Hfield_magnitude = sqrt(abs(HfieldAtPointSpherical(:, 1)).^2 + ...
%                     abs(HfieldAtPointSpherical(:,2)).^2 + ...
%                     abs(HfieldAtPointSpherical(:, 3)).^2);
%                 Hfield_magnitudedB = 20*log10(Hfield_magnitude);
%         
%                 figure;
%                 plot(r_obs, Hfield_magnitudedB);
%                 hold on
%                 xlabel('Z/m');
%                 ylabel('H/dBA/m');
        
        relError = calculateErrorNormPercentage(xVectors.Isol(1:Solver_setup.num_metallic_edges,index), Solution.PO.Isol(:,1));
        message_fc(Const,sprintf('Rel. error norm. compared to reference sol. %f percent', relError));
        
        % Calculate now the magnitude of the E-field vector.
        Efield_magnitude = sqrt(abs(EfieldAtPointSpherical(1))^2 + ...
            abs(EfieldAtPointSpherical(2))^2 + ...
            abs(EfieldAtPointSpherical(3))^2);
        RCS(index) = 4*pi*(r.*(Efield_magnitude))^2;
    end%for
end%for

end

function r = getBoundingRadius(Solver_setup)
vertices = double(Solver_setup.nodes_xyz);
r = max(sqrt(sum(vertices.^2,2)));
end











