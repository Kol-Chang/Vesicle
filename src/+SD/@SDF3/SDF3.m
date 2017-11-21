classdef SDF3 < handle
	
	%SDF3 : signed distance function in 3D

	properties (SetAccess = immutable)
		GD3 % SD.GD3 object
	end

	properties 
		F % values of the signed distance function
		Fx % 1st derivative along x-axis with central difference
		Fy % 1st derivative along y-axis with central difference
		Fz % 1st derivative along z-axis with central difference
		Fxx % 2nd derivative along x-axis with central difference
		Fyy % 2nd derivative along y-axis with central difference
		Fzz % 2nd derivative along z-axis with central difference
		Fxy % cross derivatives
		Fyz % cross derivatives
		Fzx % cross derivatives

		epsilon = 1e-4; % see 2003_Smereka
		Fg % magnitude of gradient with soomthing
		Fg_1 % magnitude of gradient without soomthing
		Nx % x component of normal vectors
		Ny % y component of normal vectors
		Nz % z component of normal vectors

		LSL % linear part of the surface Laplacian operator
		LCF % linear part of the curvature force operator
		
		SC % sum of principal curvatures. unit sphere has TC equal to 2
		CF % curvature force
		NCF % nonlinear part of the curvature force 

		LND % normal derivative operator
		CF_Op % curvature force operator
	end


	methods

		function obj = SDF3(Xm, Ym, Zm, Val)
			obj.GD3 = SD.GD3(Xm, Ym, Zm);
			obj.F = Val;
		end

		function obj = set.F(obj, val)
			obj.F = val;

			obj.Fx = (val(obj.GD3.oXo) - val(obj.GD3.oxo)) / (2*obj.GD3.Dx);
			obj.Fy = (val(obj.GD3.Yoo) - val(obj.GD3.yoo)) / (2*obj.GD3.Dy);
			obj.Fz = (val(obj.GD3.ooZ) - val(obj.GD3.ooz)) / (2*obj.GD3.Dz);
			obj.Fxx = (val(obj.GD3.oXo) - 2*val + val(obj.GD3.oxo)) / (obj.GD3.Dx.^2);
			obj.Fyy = (val(obj.GD3.Yoo) - 2*val + val(obj.GD3.yoo)) / (obj.GD3.Dy.^2);
			obj.Fzz = (val(obj.GD3.ooZ) - 2*val + val(obj.GD3.ooz)) / (obj.GD3.Dz.^2);
			obj.Fxy = (val(obj.GD3.YXo) + val(obj.GD3.yxo) - val(obj.GD3.Yxo) - val(obj.GD3.yXo)) / (4*obj.GD3.Ds.^2);
			obj.Fyz = (val(obj.GD3.YoZ) + val(obj.GD3.yoz) - val(obj.GD3.Yoz) - val(obj.GD3.yoZ)) / (4*obj.GD3.Ds.^2);
			obj.Fzx = (val(obj.GD3.oXZ) + val(obj.GD3.oxz) - val(obj.GD3.oXz) - val(obj.GD3.oxZ)) / (4*obj.GD3.Ds.^2);
			obj.Fg = sqrt(obj.Fx.^2 + obj.Fy.^2 + obj.Fz.^2 + obj.epsilon);
			obj.Fg_1 = sqrt(obj.Fx.^2 + obj.Fy.^2 + obj.Fz.^2 );
			obj.Nx = obj.Fx ./ obj.Fg;
			obj.Ny = obj.Fy ./ obj.Fg;
			obj.Nz = obj.Fz ./ obj.Fg;

			obj.LSL = obj.GD3.Lxx + obj.GD3.Lyy + obj.GD3.Lzz ...
					- ( obj.GD3.SparseDiag(obj.Nx .* obj.Nx) * obj.GD3.Lxx + ...
						obj.GD3.SparseDiag(obj.Ny .* obj.Ny) * obj.GD3.Lyy + ...
						obj.GD3.SparseDiag(obj.Nz .* obj.Nz) * obj.GD3.Lzz + ...
						obj.GD3.SparseDiag(obj.Nx .* obj.Ny) * obj.GD3.Lxy * 2 + ...
						obj.GD3.SparseDiag(obj.Ny .* obj.Nz) * obj.GD3.Lyz * 2 + ...
						obj.GD3.SparseDiag(obj.Nz .* obj.Nx) * obj.GD3.Lzx * 2  );

			obj.LND = obj.GD3.SparseDiag(obj.Nx) * obj.GD3.Lx + ...
					  obj.GD3.SparseDiag(obj.Ny) * obj.GD3.Ly + ...
					  obj.GD3.SparseDiag(obj.Nz) * obj.GD3.Lz ;


			%obj.LCF = obj.GD3.SparseDiag(obj.Fg_1) * obj.LSL * obj.GD3.SparseDiag(1./obj.Fg) * obj.LSL;
						
			obj.SC = (obj.Fxx + obj.Fyy + obj.Fzz ...
				- obj.Nx .* obj.Nx .* obj.Fxx ... 
				- obj.Ny .* obj.Ny .* obj.Fyy ...
				- obj.Nz .* obj.Nz .* obj.Fzz ...
				- obj.Nx .* obj.Ny .* obj.Fxy * 2 ...
				- obj.Ny .* obj.Nz .* obj.Fyz * 2 ...
				- obj.Nz .* obj.Nx .* obj.Fzx * 2 ) ./ obj.Fg;

			%obj.NCF = obj.Fg_1 .* obj.SC .* obj.ND(obj.SC);
			%obj.CF = obj.SL(obj.SC);

			
			
			%obj.TotalCurvature = reshape(obj.SparseDiag(1./obj.Fg) * obj.LSL * obj.F(:), obj.GD3.Size);

			obj.CF_Op = obj.GD3.SparseDiag(obj.Fg_1) * ...
					  ( obj.LSL - obj.GD3.SparseDiag(obj.SC) * obj.LND ) * ...
					  obj.GD3.SparseDiag(1./obj.Fg) * obj.LSL ;
		end

	end

	methods
		% derivative in the normal direction
		function val = ND(obj, Field)
			val = obj.Nx .* obj.GD3.Fx(Field) + obj.Ny .* obj.GD3.Fy(Field) + obj.Nz .* obj.GD3.Fz(Field);
		end

		% surface Laplacian of a Field
		function val = SL(obj, Field)
			val = obj.GD3.Laplacian(Field) - obj.SC .* obj.ND(Field) ...
				- ( obj.Nx .* obj.Nx .* obj.GD3.Fxx(Field) + ...
					obj.Ny .* obj.Ny .* obj.GD3.Fyy(Field) + ...
					obj.Nz .* obj.Nz .* obj.GD3.Fzz(Field) + ...
					obj.Nx .* obj.Ny .* obj.GD3.Fxy(Field) * 2 + ...
					obj.Ny .* obj.Nz .* obj.GD3.Fyz(Field) * 2 + ...
					obj.Nz .* obj.Nx .* obj.GD3.Fzx(Field) * 2 ) ;
		end

		function val = CurvatureForce(obj)
			val = obj.SL(obj.SC);
		end

		function val = Interp(obj, Field)
			val = griddedInterpolant(obj.GD3.PX,obj.GD3.PY,obj.GD3.PZ,permute(Field,[2,1,3]),'linear' );
		end

		function val = Surf(obj, isovalue)
			val = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F, isovalue);
		end

	end

	methods
		reinitialization(obj, Distance)
		reinitialization2(obj, Distance)
	end





	methods % visualization methods

		% plot a 3D field on the val contour of the distance function
		function plotField(obj,val,Field)
			% calculate an interpolant of the Field
			Interp = griddedInterpolant(obj.GD3.PX,obj.GD3.PY,obj.GD3.PZ,permute(Field,[2,1,3]),'linear' );
			% triangle mesh of the val isosurface. TriMesh is a structure with fields "vertices" and "faces"
			TriMesh = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F,val);
			% interpolate the values of the Field onto the vertices of the triangle mesh
			SurfField = Interp(TriMesh.vertices(:,1),TriMesh.vertices(:,2),TriMesh.vertices(:,3));
			% plot surface mesh 
			patch('Vertices',TriMesh.vertices,'Faces',TriMesh.faces,'FaceVertexCData',SurfField,...
				'FaceColor','interp','EdgeColor','none')
			axis equal
			view(45,30)
			colorbar
		end

		% plot several half contours of the distance function
		function plot(obj)
			axis(obj.GD3.BOX)
			obj.plotIso(-12*obj.GD3.Dx,0.8,'Red')
			obj.plotIso(-6*obj.GD3.Dx,0.8,'Green')
			obj.plotIso(0,0.8,'Blue')
			obj.plotIso(6*obj.GD3.Dx,0.8,'Green')
			obj.plotIso(12*obj.GD3.Dx,0.8,'Red')
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end

		% plot half of the val contour of the distance function
		function plotIso(obj,val,trans,Color)
			F = obj.F;
			F(obj.GD3.Y<0) = inf;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','none','FaceAlpha',trans);
		end

		% plot the val contour of the distance function
		function plotSurface(obj,val,trans,Color, time)
			F = obj.F;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','none','FaceAlpha',trans);
			axis(obj.GD3.BOX)
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end
	end

end