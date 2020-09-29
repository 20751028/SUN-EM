function visibility_matrix = selfShadow(Solver_setup)
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
for orig = Solver_setup.rwg_basis_functions_shared_edge_centre
    %set up rays from BF to each centeroid
    rays.orig = repmat(orig,[Solver_setup.num_metallic_triangles 1]);
    rays.dir = Solver_setup.triangle_centre_point .- orig;
    
    %Run RT to find line-of-sight visibility
    [tri_id,dist] = Raytracer.intersect_rays(rays);
    [tp, ~, itp] = unique(Solver_setup.rwg_basis_functions_trianglePlus);
    [tm, ~, itm] = unique(Solver_setup.rwg_basis_functions_triangleMinus);
    [~, i_vis_pos, ~] = intersect(tp, tri_id);
    [~, i_vis_neg, ~] = intersect(tm, tri_id);
    tp(i_vis_pos) = -1;
    vis_pos = (tp(itp)==-1);
    tp(i_vis_neg) = -1;
    vis_neg = (tp(itm)==-1);
    visible = vis_pos & vis_neg;
    visible = visible(1:Solver_setup.num_metallic_triangles);
end


end