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
r = 1000;
% Calculate now the E-field value here internal
index = 0;
RCS = zeros(length(theta_grid)*length(phi_grid), 1);
for theta_degrees = theta_grid
    for phi_degrees = phi_grid % 0 to 90 degr. in steps of 1 degrees
        index = index + 1;
        
        Solver_setup.theta = theta_degrees;
        Solver_setup.phi = phi_degrees;
        
        % --------------------------------------------------------------------------------------------------
        % Run the EM solver
        % --------------------------------------------------------------------------------------------------
        [Solution] = runEMsolvers(Const, Solver_setup, 0, 0, xVectors);
        
        EfieldAtPointSpherical =  calculateEfieldAtPointRWG(Const, r, theta_degrees, phi_degrees, ...
            Solver_setup, Solution.PO.Isol);
        
        
        % Calculate now the magnitude of the E-field vector.
        Efield_magnitude = sqrt(abs(EfieldAtPointSpherical(1))^2 + ...
            abs(EfieldAtPointSpherical(2))^2 + ...
            abs(EfieldAtPointSpherical(3))^2);
        RCS(index) = 4*pi*(r.*(Efield_magnitude))^2;
    end%for
end%for

end











