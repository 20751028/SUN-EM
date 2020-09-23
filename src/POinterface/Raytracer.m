classdef Raytracer
    methods(Static)
        function setGeom(Solver_setup)
            
            face_index = Solver_setup.triangle_vertices-1;
            
            % Send geometry to raytracer
            % Note: mexRaySrv creates and preserves a copy of the vertex
            %       and face index data, so we don't have to preserve the
            %       variables here.
            vertex_coord = Solver_setup.nodes_xyz;
            modified = datestr(now);
            mexRaySrv('setmesh',vertex_coord,face_index, modified);
            
        end
        function [tri_id,dist] = intersect_rays(rays)
            % Find ray intersections
            mexRaySrv('castrays',[rays.origin, rays.dir]);
            [tri_id,dist] = mexRaySrv('getresults');
            tri_id = tri_id+1;
        end
    end
end