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
Const.FEKOmatfilename          = 'plate2.mat';
Const.FEKOstrfilename          = 'plate2.str';
Const.FEKOrhsfilename          = 'plate2.rhs';
Const.FEKOoutfilename          = 'plate2.out';
Const.FEKOefefilename          = 'plate2_NearField1.efe';
Const.FEKOffefilename          = 'plate2.ffe';

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
theta_grid = 0:1:90;
phi_grid = 45:1:45;

num_theta_samples = length(theta_grid);
num_phi_samples = length(phi_grid);
total_efield_samples = num_theta_samples*num_phi_samples;
Efield_magnitude = zeros(total_efield_samples,1);

%calculate the monostatic RCS
RCS = calcRCS(Const, Solver_setup, theta_grid, phi_grid, xVectors);

% Plot now the RCS
figure;
hold on;
grid on;

plot(1:total_efield_samples,10*log10(RCS));
set(get(gca, 'XLabel'), 'String', ('Sample index'));
set(get(gca, 'YLabel'), 'String', ('RCS [dBsm]'));

