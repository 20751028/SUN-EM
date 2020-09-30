function [po] = runPOsolver(Const, Solver_setup, ~, ~, refIsol)

%runPOsolver
%   Usage:
%       [po] = runPOsolver(Const)
%
%   Input Arguments:
%       Const
%           A global struct, containing general data
%       Solver_setup
%           Solver specific struct, e.g. frequency range, basis function details, geometry details
%       zMatrices
%           The Z-matrices data. This can be from FEKO (extracted from the *.mat file, or internally
%           calculated).
%       yVectors
%           The Yrhs-vector data
%       refIsol
%           The reference solution-vector data (e.g. MoM solution of FEKO or SUN-EM)
%
%   Output Arguments:
%       po
%           Structs containing PO solution and timing data
%
%   Description:
%       Runs the PO solution based on the Z and Y data that was read / parsed
%       from the FEKO *.out, *.mat, *.str and *.rhs files
%
%   =======================
%   Written by Cullen Stewart-Burger on 27 July 2020
%   Adapted from runMoMsolver.m by Danie Ludick
%   Stellenbosch University
%   Email: 20751028@sun.ac.za

%   indev notes (28/07/2020):
%       Current implementation assumes fully illuminated target
%       Excitation source is set to a unit phi-polarised plane wave
%       travelling in -r direction.

narginchk(5,5);

message_fc(Const,' ');
message_fc(Const,'------------------------------------------------------------------------------------');
message_fc(Const,sprintf('Running PO solver'));

% Initialise the return values
po  = [];
po.name = 'po';
Npo = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
%numSols = xVectors.numSols;                   % The number of solutions configurations
po.numSols = 1; %numSols;                     % For now, set to 1. (TO-DO: Update)
numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
numRHSperFreq = po.numSols / numFreq;         % The number of solutions per frequency point.
% For now, should be 1 (TO-DO: Update)

% Some info about the solution configurations
% message_fc(Const,sprintf('  numSols : %d', po.numSols));
% message_fc(Const,sprintf('  numFreq : %d', numFreq));
% message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));

% Calculate the solution vector (observable BFs only)
po.Isol = complex(zeros(Npo,1));
%create an intermediate solution vector with all BFs positive and negative
Isol = complex(zeros(Npo,2));
% The timing calculations also need to take into account that there is a
% frequency loop
po.setupTime = zeros(1,numFreq);
% Zero also the total times (for all frequency iterations)
po.totsetupTime = 0.0;
po.totfactorisationTime = 0.0;
po.totsolTime = 0.0;

if(Solver_setup.num_reflections > 1)
    MRPO = true;
end


%Test to see which basis functions are illuminated
ray = setRays(Solver_setup);
[tri_id, dist] = Raytracer.intersect_rays(ray);
[tp, ~, itp] = unique(Solver_setup.rwg_basis_functions_trianglePlus);
[tm, ~, itm] = unique(Solver_setup.rwg_basis_functions_triangleMinus);
[~, i_vis_pos, ~] = intersect(tp, tri_id);
[~, i_vis_neg, ~] = intersect(tm, tri_id);
tp(i_vis_pos) = -1;
vis_pos = (tp(itp)==-1);
tp(i_vis_neg) = -1;
vis_neg = (tp(itm)==-1);
visible = vis_pos & vis_neg;
visible = visible(1:Npo);
%debug: plot the visible basis functions
%scatter3(Solver_setup.rwg_basis_functions_shared_edge_centre(visible,1),Solver_setup.rwg_basis_functions_shared_edge_centre(visible,2),Solver_setup.rwg_basis_functions_shared_edge_centre(visible,3));
%axis('equal');

% Start the frequency loop now
for freq=1:numFreq
    
    % Start timing (PO setup)
    tic
    
    % Reset each solution per frequency point here
    solStart = 1;
    solEnd   = numRHSperFreq;
    
    %         % Allocate here space for the MoM matrix so that it can be filled
    %         ObservRWGs = [1:Npo];
    %         SourceRWGs = [1:Npo];
    %         % Note: Since 2017-06-25, we are also passing a freq. variable here to
    %         % indicate for which frequency point we are extracting the matrix
    %         Zmom = (calcZmn(Const, zMatrices, freq, 1, 1, ObservRWGs, SourceRWGs));
    
    % End timing (calculating the impedance matrix)
    po.setupTime(freq) = toc;
    
    % Start timing
    tic
    
    
    rn = Solver_setup.rwg_basis_functions_shared_edge_centre;
    shared_nodes = Solver_setup.rwg_basis_functions_shared_edge_nodes;
    ln = Solver_setup.nodes_xyz(shared_nodes(:, 2), :) - Solver_setup.nodes_xyz(shared_nodes(:, 1), :);
    ln = ln./Solver_setup.rwg_basis_functions_length_m(1:Npo);
    %We do not know if this direction for ln is correct according
    %to our reference. This is checked and corrected if necessary
    %below.
    
    side = Solver_setup.nodes_xyz(Solver_setup.rwg_basis_functions_trianglePlusFreeVertex, :) - Solver_setup.nodes_xyz(shared_nodes(:, 1), :);
    normTest = cross(side, ln, 2);
    reverse = dot(normTest, Solver_setup.triangle_normal_vector(Solver_setup.rwg_basis_functions_trianglePlus(1:Npo), :), 2);
    ln = ln.*reverse./abs(reverse);
    
    
    %The visibility term (delta) should be set to 0 if either triangle is not
    %visible or +-1 depending on the incident direction relative to the
    %surface normal.
    nDotRay = dot(Solver_setup.triangle_normal_vector(Solver_setup.rwg_basis_functions_trianglePlus(1:Npo), :), ray.dir(Solver_setup.rwg_basis_functions_trianglePlus(1:Npo), :), 2);
    delta = nDotRay./abs(nDotRay).*visible;
    BF_side = [delta > 0, delta < 0];
    
    k = 2*pi*Solver_setup.frequencies.samples(freq)/Const.C0;
    [kx, ky, kz] = sph2cart(Solver_setup.phi*Const.DEG2RAD, (90-Solver_setup.theta)*Const.DEG2RAD, -k);
    k_vec = repmat([kx, ky, kz], [Npo, 1]);
    H = (1/Const.ETA_0)*exp(j*dot(k_vec, rn, 2));    %Find impressed H field (r directed plane wave with E field theta-polarised)
    a_phi = [-sind(Solver_setup.phi), cosd(Solver_setup.phi), 0];
    H_vec = -H*a_phi;
    Isol = repmat(dot(2*delta.*H_vec, ln, 2), [1, 2]).*BF_side;
    
    
    for refl_num = 2:Solver_setup.num_reflections
        Isol_refl = complex(zeros(Npo, 2));
        tic
        for n = 1:Npo
            Isol_pos = zeros(Npo, 1);
            Isol_neg = zeros(Npo, 1);
            Isol_pos((Solver_setup.Visibility_matrix(:, n) == 1)) = Isol((Solver_setup.Visibility_matrix(:, n) == 1), 1);
            Isol_pos((Solver_setup.Visibility_matrix(:, n) == 3)) = Isol((Solver_setup.Visibility_matrix(:, n) == 3), 2);
            Isol_neg((Solver_setup.Visibility_matrix(:, n) == 2)) = Isol((Solver_setup.Visibility_matrix(:, n) == 2), 1);
            Isol_neg((Solver_setup.Visibility_matrix(:, n) == 4)) = Isol((Solver_setup.Visibility_matrix(:, n) == 4), 2);
            H_vec_pos = calculateHfieldAtPointRWGCart(Const, rn(n, :), Solver_setup, Isol_pos);
            %H_vec_neg = calculateHfieldAtPointRWG(Const, rn(n, 1), rn(n, 3), rn(n, 3), Solver_setup, Isol_neg);
            H_vec_neg = calculateHfieldAtPointRWGCart(Const, rn(n, :), Solver_setup, Isol_neg);
            Isol_refl(n, :) = [dot(2*H_vec_pos, ln(n, :), 2)  dot(-2*H_vec_neg, ln(n, :), 2)];
        end
        toc
        Isol = Isol + Isol_refl;
    end
    
    % End timing (MoM factorisation)
    po.factorisationTime(freq) = toc;
    
    % Total time (MoM matrix setup + factorisation)
    po.solTime(freq) = po.setupTime(freq) + po.factorisationTime(freq);
    
    % Memory usage remains constant between frequency iterations
    po.memUsage = 0;%byteSize(Zmom);
    
    % Calculate the total MoM times
    po.totsetupTime = po.totsetupTime + po.setupTime(freq);
    po.totfactorisationTime = po.totfactorisationTime + po.factorisationTime(freq);
    po.totsolTime = po.totsolTime + po.solTime(freq);
    
end%for freq=1:numFreq
po.Isol = po.Isol + Isol(:, 1).*BF_side(:, 1);
po.Isol = po.Isol + Isol(:, 2).*BF_side(:, 2);
%po.Isol = Isol(:, 1) + Isol(:, 2);

message_fc(Const,sprintf('Finished PO solver in %f sec.',po.totsolTime));

% Compare the MoM solution obtained with MATLAB, with that obtained by FEKO
% that was stored in xVectors.values (for each frequency iteration (and each solution within the frequency iteration)
% Calculate also space for the relative error here
% po.relError = zeros(1,po.numSols);
% for freq=1:numFreq
%     for solNum=1:numRHSperFreq
%         index = solNum + (freq-1)*numRHSperFreq;
%         po.relError(index) = calculateErrorNormPercentage(refIsol.Isol(1:Solver_setup.num_metallic_edges,index), po.Isol(:,index));
%         message_fc(Const,sprintf('Rel. error norm. for Sol. %d of %d of freq. %d of %d compared to reference sol. %f percent',solNum, ...
%             numRHSperFreq, freq, numFreq, po.relError(index)));
%     end
% end

% Write the MoM solution to a ASCII str file, so that it can be read
% again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
if (~isempty(Const.SUNEMmomstrfilename))
    writeSolToFile(Const, po);
end%if

