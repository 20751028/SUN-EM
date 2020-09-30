function Visibility_matrix = selfShadow(Solver_setup)
% for closed surfaces:
% for each bf
%     for each centroid
%         set up ray from bf to centroid
%     cast all rays for given bf
%     check if ray intersects with corresp. triangle, tri is visible
%     dot ray with the triange normal to check if the ray passed through the solid object (can this be done before ray casting?)
%     test if pos and neg triangles are visible
%     bf is visible if both triangles are visible
%     
%  for generalised surfaces, we need something more sophisticated
%      consider replicating all the geometry as "inner" and "outer" surfaces
%      in this case, we need to follow the process above (except RT which can be done once) with four seperate cases:
%      !!check below!!
%      outside BF to outside centroids - ray and normal same direction angle > 90 = solid
%      outside BF to inside centroid - ray and normal same dir = solid
%      inside BF to outside centroid - ray and normal 
%      inside BF to inside centroid - ray and normal

%set up visibility matrix functions
% Vpp = logical(sparse(Solver_setup.num_mom_basis_functions, Solver_setup.num_mom_basis_functions));
% Vpn = logical(sparse(Solver_setup.num_mom_basis_functions, Solver_setup.num_mom_basis_functions));
% Vnp = logical(sparse(Solver_setup.num_mom_basis_functions, Solver_setup.num_mom_basis_functions));
% Vnn = logical(sparse(Solver_setup.num_mom_basis_functions, Solver_setup.num_mom_basis_functions));
Visibility_matrix = zeros(Solver_setup.num_mom_basis_functions, Solver_setup.num_mom_basis_functions, 'uint8');

%Before running visibility tests, we need to construct an outward facing
%normal vector for each BF, as defined in Xiang, 2016.
rwg_normal_vec = zeros(Solver_setup.num_mom_basis_functions, 3);
for ind = 1:Solver_setup.num_mom_basis_functions
rwg_normal_vec(ind, :) = Solver_setup.triangle_normal_vector(Solver_setup.rwg_basis_functions_trianglePlus(ind), :)...
    + Solver_setup.triangle_normal_vector(Solver_setup.rwg_basis_functions_triangleMinus(ind), :);
end
%note: We have not normalised the vector, but as we are only interested in
%direction this is not necesssary.

for ind = 1:Solver_setup.num_mom_basis_functions
    orig = Solver_setup.rwg_basis_functions_shared_edge_centre(ind, :);
    %set up rays from BF to each centeroid
    rays.origin = repmat(orig,[Solver_setup.num_metallic_triangles 1]);
    rays.dir = Solver_setup.triangle_centre_point - rays.origin;
    
    %Run RT to find line-of-sight visibility
    [tri_id,~] = Raytracer.intersect_rays(rays);
    [tp, ~, itp] = unique(Solver_setup.rwg_basis_functions_trianglePlus);
    [tm, ~, itm] = unique(Solver_setup.rwg_basis_functions_triangleMinus);
    [~, i_vis_pos, ~] = intersect(tp, tri_id);
    [~, i_vis_neg, ~] = intersect(tm, tri_id);
    tp(i_vis_pos) = -1;
    vis_pos = (tp(itp)==-1);
    tp(i_vis_neg) = -1;
    vis_neg = (tp(itm)==-1);
    visible = vis_pos & vis_neg;
    visible = visible(1:Solver_setup.num_mom_basis_functions);
    %we now have a list of all the basis functions with line-of-sight
    %visibility from the source basis function centre point
    
    %Next, we need to determine which sides of the BFs are facing each other
    r_source_obs = Solver_setup.rwg_basis_functions_shared_edge_centre -...
        repmat(orig,[Solver_setup.num_mom_basis_functions 1]);
    source_norm = repmat(rwg_normal_vec(ind, :), [Solver_setup.num_mom_basis_functions 1]);
    source_facing = dot(source_norm, r_source_obs, 2);
    obs_facing = dot(rwg_normal_vec, r_source_obs, 2);
    
    pp = (source_facing > 0) & (obs_facing < 0);
    pn = (source_facing > 0) & (obs_facing > 0);
    np = (source_facing < 0) & (obs_facing < 0);
    nn = (source_facing < 0) & (obs_facing > 0);
    
    Visibility_matrix(ind, pp) = visible(pp)*1;
    Visibility_matrix(ind, pn) = visible(pn)*2;
    Visibility_matrix(ind, np) = visible(np)*3;
    Visibility_matrix(ind, nn) = visible(nn)*4;
end
end