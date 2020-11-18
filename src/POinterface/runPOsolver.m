function [po] = runPOsolver(Const, Solver_setup, ~, ~, refIsol)

%runPOsolver
%   Usage:
%       [po] = runPOsolver(Const)
%
%   Input Arguments:
%       Const
%           A global struct, containing general data
%       Solver_setup
%           Solver specific struct, e.g. frequency range, basis function
%           details, geometry details, aspect angles
%       refIsol
%           The reference solution-vector data (e.g. PO solution of FEKO or SUN-EM)
%
%   Output Arguments:
%       po
%           Structs containing PO solution and timing data
%
%   Description:
%       Calculates the PO solution based on the Z and Y data that was read / parsed
%       from the FEKO *.out fiel
%
%   =======================
%   Written by Cullen Stewart-Burger on 27 July 2020
%   Stellenbosch University
%   Email: 20751028@sun.ac.za

narginchk(5,5);

message_fc(Const,' ');
message_fc(Const,'------------------------------------------------------------------------------------');
message_fc(Const,sprintf('Running PO solver'));

% Initialise the return values
po  = [];
po.name = 'po';
Npo = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
po.numSols = 1; %numSols;                     % For now, set to 1. (TO-DO: Update)
numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
numRHSperFreq = po.numSols / numFreq;         % The number of solutions per frequency point.

% Calculate the solution vector (observable BFs only)
po.Isol = complex(zeros(Npo,1));
%create an intermediate solution vector with all BFs positive and negative
Isol = complex(zeros(Npo,2));
% The timing calculations also need to take into account that there is a
% frequency loop
po.setupTime = zeros(1,numFreq);
% Zero also the total times (for all frequency iterations)
po.totsetupTime = 0.0;
po.totsolTime = 0.0;


%Test to see which basis functions are illuminated
ray = setRays(Solver_setup, Solver_setup.theta, Solver_setup.phi);
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
Isol = complex(zeros(numFreq, Npo));
for freq=1:numFreq
    
    % Start timing (PO setup)
    tic
    
    % Reset each solution per frequency point here
    solStart = 1;
    solEnd   = numRHSperFreq;
    
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
    delta = -nDotRay./abs(nDotRay).*visible;
    %delta = -nDotRay./abs(nDotRay); %uncomment this line to use the full
    %illumination assumption
    BF_side = [delta > 0.017, delta < -0.017];
    delta(isnan(delta)) = 0;
    
    k = 2*pi*Solver_setup.frequencies.samples(freq)/Const.C0;
    [kx, ky, kz] = sph2cart(Solver_setup.phi*Const.DEG2RAD, (90-Solver_setup.theta)*Const.DEG2RAD, -k);
    k_vec = repmat([kx, ky, kz], [Npo, 1]);
    H = (1/Const.ETA_0)*exp(1j*(dot(k_vec, rn, 2)));    %Find impressed H field (-r directed plane wave with E field theta-polarised)
    a_phi = [-sind(Solver_setup.phi), cosd(Solver_setup.phi), 0];
    a_theta = -[cosd(Solver_setup.theta)*cosd(Solver_setup.phi), cosd(Solver_setup.theta)*sind(Solver_setup.phi), -sind(Solver_setup.theta)];
    H_vec = H*a_phi;
    Isol = repmat(dot(2*delta.*H_vec, ln, 2), [1, 2]).*BF_side;
    %Isol(freq, :) = dot(2*delta.*H_vec, ln, 2);
    
    Isol_refl = complex(zeros(Npo, 2, Solver_setup.num_reflections));
    Isol_refl(:, :, 1) = Isol;
    for refl_num = 2:Solver_setup.num_reflections
        H_vec_pos = zeros(Npo, 3);
        H_vec_neg = zeros(Npo, 3);
        for n = 1:Npo
            Isol_pos = complex(zeros(Npo, 1));
            Isol_neg = complex(zeros(Npo, 1));
            Isol_pos((Solver_setup.Visibility_matrix(:, n) == 1)) = Isol_refl((Solver_setup.Visibility_matrix(:, n) == 1), 1, refl_num-1);
            Isol_pos((Solver_setup.Visibility_matrix(:, n) == 3)) = Isol_refl((Solver_setup.Visibility_matrix(:, n) == 3), 2, refl_num-1);
            Isol_neg((Solver_setup.Visibility_matrix(:, n) == 2)) = Isol_refl((Solver_setup.Visibility_matrix(:, n) == 2), 1, refl_num-1);
            Isol_neg((Solver_setup.Visibility_matrix(:, n) == 4)) = Isol_refl((Solver_setup.Visibility_matrix(:, n) == 4), 2, refl_num-1);
            H_vec_pos(n, :) = calculateHfieldAtPointRWGCart(Const, rn(n, :), Solver_setup, Isol_pos);
            H_vec_neg(n, :) = calculateHfieldAtPointRWGCart(Const, rn(n, :), Solver_setup, Isol_neg);
        end
        Isol_refl(:, :, refl_num) = [dot(2*conj(H_vec_pos), ln, 2)  dot(-2*conj(H_vec_neg), ln, 2)];
    end
    Isol = sum(Isol_refl, 3);
    
    
    % Memory usage remains constant between frequency iterations
    %po.memUsage = byteSize(Solver_setup.Visibility_matrix);
    
    
end%for freq=1:numFreq

if(Solver_setup.is_bistatic)
    %Test to see which basis functions are visible from the reciever
    theta = Solver_setup.theta + Solver_setup.theta_bistatic;
    phi = Solver_setup.phi + Solver_setup.phi_bistatic;
    ray = setRays(Solver_setup, theta, phi);
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
    
    %Test which side of each basis function is visible
    nDotRay = dot(Solver_setup.triangle_normal_vector(Solver_setup.rwg_basis_functions_trianglePlus(1:Npo), :), ray.dir(Solver_setup.rwg_basis_functions_trianglePlus(1:Npo), :), 2);
    delta = -nDotRay./abs(nDotRay).*visible;
    BF_side = [delta > 0.017, delta < -0.017];
end


po.Isol = Isol(:, 1).*BF_side(:, 1);
po.Isol = po.Isol + Isol(:, 2).*BF_side(:, 2);
po.totsolTime = toc;

message_fc(Const,sprintf('Finished PO solver in %f sec.', po.totsolTime));



