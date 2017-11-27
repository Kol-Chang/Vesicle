Kappa = 411;

% vesicle dynamics
[x, y, z] = meshgrid(linspace(-250,250,64)); % simulation domain in nm

a = 215;
b = 215;
c = 70;
F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;
%F = sqrt(x.^2+y.^2+z.^2) - 215;
map = SD.SDF3(x,y,z,F);
map.reinitialization( map.F )
map.plotSurface(0,1,'g')

loops = 100;
Skip = 1;
SkipR = 1;
Dt = 1 * map.GD3.Dx ^ 4;

for ii = 1:loops-1

	% calculate normal velocity due to bending
	co = 0.5 * (map.SC.^2 - 4*map.GC); % will be reused
	Vb = Kappa * (map.SL(map.SC) + map.SC .* co);
	% change of area and volume due to Vb
	DA = map.SurfaceIntegral(map.SC .* Vb);
	DV = map.SurfaceIntegral(Vb);
	s_mat = [DA; DV];
	% coefficient matrix from lagrange multiplier
	c11 = map.SurfaceIntegral(map.SC.^2);
	c12 = - map.SurfaceIntegral(map.SC);
	c21 = - c12;
	c22 = - map.Suf_Area;
	c_mat = [c11, c12; c21, c22]; 
	% calculate surface tension and pressure
	lag_mul = c_mat\s_mat; % lagrange multiplier
	Tension = lag_mul(1);
	Pressure = lag_mul(2);

	%% construct the force density operator
	% an operator appearing in surface laplacian operator and curvature operator
	LOP = map.GD3.Lxx + map.GD3.Lyy + map.GD3.Lzz ...
				- ( map.GD3.SparseDiag(map.Nx .* map.Nx) * map.GD3.Lxx + ...
					map.GD3.SparseDiag(map.Ny .* map.Ny) * map.GD3.Lyy + ...
					map.GD3.SparseDiag(map.Nz .* map.Nz) * map.GD3.Lzz + ...
					map.GD3.SparseDiag(map.Nx .* map.Ny) * map.GD3.Lxy * 2 + ...
					map.GD3.SparseDiag(map.Ny .* map.Nz) * map.GD3.Lyz * 2 + ...
					map.GD3.SparseDiag(map.Nz .* map.Nx) * map.GD3.Lzx * 2  ); 
	% curvature operator 	
	CvL = map.GD3.SparseDiag(1./map.Fg) * LOP;
	% normal derivative
	ND = map.GD3.SparseDiag(map.Nx) * map.GD3.Lx + ...
		 map.GD3.SparseDiag(map.Ny) * map.GD3.Ly + ...
		 map.GD3.SparseDiag(map.Nz) * map.GD3.Lz ;
	% surface laplacian operator
	SL = LOP - map.GD3.SparseDiag(map.SC) * ND;
	% normal velocity operator without drag coefficient and Pressure
	%NV = Kappa * ( ( SL + map.GD3.SparseDiag(co) ) * CvL ) - Tension * CvL;
	NV = Kappa * ( ( SL + map.GD3.SparseDiag(co) ) * CvL );

	% the advection term without pressure
	A = map.GD3.SparseDiag(map.Fg_1) * NV * Dt / 2;
	B = map.GD3.Idt + A;
	C = map.GD3.Idt - A;

	F_old = map.F;
	%S = C * F_old(:) - Dt * Pressure;
	S = C * F_old(:);

	[L,U]=ilu(B,struct('type','nofill','milu','row'));
	F_new = gmres(B, S, 50, 1e-12, 10, L, U);
	map.F = reshape(F_new, map.GD3.Size);

	clf
	map.plotSurface(0,1,'g')
	time = num2str(ii*Dt);
	title([num2str(ii) ': ' num2str(ii*Dt)])
	text(map.GD3.xmin,map.GD3.ymax,(map.GD3.zmax+map.GD3.zmin)/2,['BR',num2str(ii),':',time])
	drawnow

	if (mod(ii,SkipR)==0)
		map.reinitialization( reshape(F_new, map.GD3.Size) );
	end




end