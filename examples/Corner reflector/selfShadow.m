function selfShadow(Solver_setup)
for closed surfaces:
for each bf
    for each centroid
        set up ray from bf to centroid
    cast all rays for given bf
    check if ray intersects with corresp. triangle, tri is visible
    dot ray with the triange normal to check if the ray passed through the solid object (can this be done before ray casting?)
    test if pos and neg triangles are visible
    bf is visible if both triangles are visible
    
 for generalised surfaces, we need something more sophisticated
     consider replicating all the geometry as "inner" and "outer" surfaces
     in this case, we need to follow the process above (except RT which can be done once) with four seperate cases:
     !!check below!!
     outside BF to outside centroids - ray and normal same direction angle > 90 = solid
     outside BF to inside centroid - ray and normal same dir = solid
     inside BF to outside centroid - ray and normal 
     inside BF to inside centroid - ray and normal