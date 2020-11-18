% Author: Cullen Stewart-Burger (20751028@sun.ac.za)
% Project: RCS prediction of a flat plate
%
% Note: Each project directory / example directory needs to have a sunem_initialise.m
% script, that is used to setup the correct environment for that example.
%
% Refer to the /doc folder for more information

% --------------------------------------------------------------------------------------------------
% Initialise the environment
% --------------------------------------------------------------------------------------------------
% Project output directory: './results/'
% Debug: True/False
Const = sunem_initialise('results',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runPOsolver              = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOstrfilename          = 'reflectorRT1G_rev.str';
Const.FEKOoutfilename          = 'reflectorRT1G_rev.out';

% The Following file is used to port solutions to FEKO 
% (for post-processing in POSTFEKO).
% TO-DO: [DL] Add this.
% Const.output_strfilename    = '';
% Const.writeFEKOstrfile = [0 0 0 0];

% --------------------------------------------------------------------------------------------------
% Read the MoM matrix equation from the file
% --------------------------------------------------------------------------------------------------
[xVectors] = readFEKOXvectorFromFile(Const, Const.FEKOstrfilename);
yVectors.numRhs = xVectors.numMoMbasis;
Const = sunem_init(Const, yVectors);

% --------------------------------------------------------------------------------------------------
% Parse the setup files to extract the frequency sweep, the geometry and basis function setup 
% --------------------------------------------------------------------------------------------------
% TO-DO: At a later stage we can also add other meshing / geometry
% preprocessxing, e.g. Gmsh or GiD. For now the solver setup is read from FEKO.
[Const, Solver_setup] = parseFEKOoutfile(Const, 0);


% --------------------------------------------------------------------------------------------------
% Set up a number of frames to be run (ie one per incident angle)
% --------------------------------------------------------------------------------------------------
theta_grid = 10:0.5:170;
phi_grid = 0:1:0;
%Setup for bistatic case
Solver_setup.is_bistatic = true;
Solver_setup.phi_bistatic = 0;
Solver_setup.theta_bistatic = -10;


%Set number of reflections for MRPO
Solver_setup.num_reflections = 3;

num_theta_samples = length(theta_grid);
num_phi_samples = length(phi_grid);
total_efield_samples = num_theta_samples*num_phi_samples;
Efield_magnitude = zeros(total_efield_samples,1);
%calculate the monostatic RCS
RCS = calcRCS(Const, Solver_setup, theta_grid, phi_grid, xVectors);
RCS_mag = sum(RCS, 2);

% Plot the RCS
figure;
hold on;
grid on;

plot(1:total_efield_samples,10*log10(RCS_mag));
set(get(gca, 'XLabel'), 'String', ('Sample index'));
set(get(gca, 'YLabel'), 'String', ('RCS [dBsm]'));

