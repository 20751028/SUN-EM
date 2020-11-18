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
%           The reference solution-vector data (e.g. PO solution of FEKO or SUN-EM)
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
RCS = zeros(length(theta_grid)*length(phi_grid), 2);
Raytracer.setGeom(Solver_setup);
if(Solver_setup.num_reflections > 1)
    Solver_setup.Visibility_matrix = selfShadow(Solver_setup); 
    memUsage = byteSize(Solver_setup.Visibility_matrix);
end
for theta_degrees = theta_grid
    for phi_degrees = phi_grid
        
        Solver_setup.theta = theta_degrees;
        Solver_setup.phi = phi_degrees;
        
        % --------------------------------------------------------------------------------------------------
        % Run the EM solver
        % --------------------------------------------------------------------------------------------------
        index = index + 1;
        [Solution] = runEMsolvers(Const, Solver_setup, 0, 0, xVectors);
        
        %Change reviever position for bistatic case
        if(Solver_setup.is_bistatic)
            reciever_theta_degrees = theta_degrees + Solver_setup.theta_bistatic;
            reciever_phi_degrees = phi_degrees + Solver_setup.phi_bistatic;
        else
            reciever_theta_degrees = theta_degrees;
            reciever_phi_degrees = phi_degrees;
        end
        
        EfieldAtPointSpherical =  calculateEfieldAtPointRWG(Const, r, reciever_theta_degrees, reciever_phi_degrees, ...
            Solver_setup, Solution.PO.Isol);
        
        relError = calculateErrorNormPercentage(xVectors.Isol(1:Solver_setup.num_metallic_edges,index), Solution.PO.Isol(:,1));
        message_fc(Const,sprintf('Rel. error norm. compared to reference sol. %f percent', relError));
        
        % Calculate now the magnitude of the E-field vector.
        Efield_magnitude = sqrt(abs(EfieldAtPointSpherical(1))^2 + ...
            abs(EfieldAtPointSpherical(2))^2 + ...
            abs(EfieldAtPointSpherical(3))^2);
        E_theta = abs(EfieldAtPointSpherical(2));
        E_phi = abs(EfieldAtPointSpherical(3));
        RCS(index, 1) = 4*pi*(r.*(E_theta))^2;
        RCS(index, 2) = 4*pi*(r.*(E_phi))^2;
       % end
    end%for
end%for

end

function r = getBoundingRadius(Solver_setup)
vertices = double(Solver_setup.nodes_xyz);
r = max(sqrt(sum(vertices.^2,2)));
end