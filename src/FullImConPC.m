% before running a test, gather necessary information about this test
% so that the result can be recovered easily

% time of excution
	Time_Now = datetime('now');
	FormatOut = 'yy_mm_dd_HH_MM_SS_';
	Time_Begin = datestr(Time_Now,FormatOut);

% system info
	if ismac
		PlatForm = 'Mac_';
	elseif isunix
		PlatForm = 'unix_';
	elseif ispc
		PlatForm = 'Windows_';
	else
		PlatForm = 'unknown_';
	end

% make a directory tagged with time of excution outsider src folder
	% path of this folder. we randomized the last 3 digits to avoid name clashes
	Result_Folder = fullfile('..','exc',[PlatForm Time_Begin num2str(floor(rand()*1000))]); 
	if ~exist(fullfile('..','exc')) % if ../exc does not exist
		mkdir(fullfile('..','exc')) % make this folder
	end
	while exist(Result_Folder) % if this folder name alread exists -- which will probably never happen
		Result_Folder = fullfile('..','exc',[PlatForm Time_Begin num2str(floor(rand()*1000))]); % rename it
	end
	mkdir(Result_Folder) % make a folder
	
% now copy the src code the Result_Folder
	Current_src = fullfile(Result_Folder,'src'); % folder to hold current src files
	mkdir(Current_src)
	copyfile(fullfile('..','src'),Current_src);

% now open a text file and write to it comments and testing info
	Test_info = fopen(fullfile(Result_Folder,'test_info'), 'w');
	fprintf(Test_info, 'test start time : \t');
	fprintf(Test_info, [datestr(Time_Now, 'yy/mm/dd HH:MM:SS'),'\n']);
	fprintf(Test_info, 'test file: FullImConPC \n');
	fprintf(Test_info, 'test info: try to accerlerate convergence of iters \n');
	fclose(Test_info);

D = 0.5;
Alpha = 2;

[x, y, z] = meshgrid(linspace(-pi-D,pi+D,64)/Alpha);

C1 = pi/(2*Alpha); 
%C2 = 0.90 * C1/2;
C2 = 0.8 * C1/2;

F1 = sqrt(x.^2+y.^2) - (C1-C2*(cos(Alpha * z)+1));
F2 = max(z-pi/Alpha,-z-pi/Alpha);
F = max(F1, F2);

map = SD.SDF3(x,y,z,F)
map.reinitialization(F)

% save grids
GridX = map.GD3.X;
GridY = map.GD3.Y;
GridZ = map.GD3.Z;
save(fullfile(Result_Folder,'Grid.mat'),'GridX','GridY','GridZ');

Dt = 20 * map.GD3.Dx ^ 4;

loops = 100;
Skip = 1;
SkipR = 1;
count= 1;
mov(ceil(loops/Skip)) = struct('cdata',[],'colormap',[]);
%snap{1000} = [];

figure

for ii = 1:loops-1

	A = map.CF_Op * Dt /2;
	B = map.GD3.Idt + A;
	C = map.GD3.Idt - A;



	F_old = map.F;
	S = C * F_old(:);

	keyboard

	%F_new = bicg(B, S);
	%F_new = bicgstab(B, S, [], 50);
	[L,U]=ilu(B,struct('type','nofill','milu','row')); % this seems to be the only pratical setup 

	% preconditioner

	% it seems that providing a intial guess will slow the algorithm down with either methods
	% after some comparing bicgstab outperforms all other mehtods except gmres wihtout restart
	% gmres with restart less than 30 performs worse than bicgstab
	% with no restart gmres provide best performance
	%F_new = bicgstab(B, S, 1e-12, 200, L, U);
	F_new = gmres(B, S, 50, 1e-12, 10, L, U);
	%F_new = pcg(B, S);
	%map.F = reshape(F_new, map.GD3.Size);
	map.F = reshape(F_new, map.GD3.Size);

	clf
	map.plotSurface(0,1,'g')
	time = num2str(ii*Dt);
	title(num2str(ii*Dt))
	text(map.GD3.xmin,map.GD3.ymax,(map.GD3.zmax+map.GD3.zmin)/2,['BR',num2str(ii),':',time])
	drawnow

	%saveas(gcf, fullfile(Result_Folder,[num2str(ii),'AA_BR','.png']))

	if (mod(ii,Skip)==0)
	%	mov(count) = getframe(gcf);
		DistanceMap = map.F;
		saveas(gcf, fullfile(Result_Folder,[num2str(ii),'AA_BR','.png']))
		save(fullfile(Result_Folder,['DFV',num2str(ii),'AA_BR','.mat']),'DistanceMap')
	%	count = count + 1;
	end

	if (mod(ii,SkipR)==0)
		map.reinitialization( reshape(F_new, map.GD3.Size) );
		clf
		map.plotSurface(0,1,'g')
		time = num2str(ii*Dt);
		title(num2str(ii*Dt))
		text(map.GD3.xmin,map.GD3.ymax,(map.GD3.zmax+map.GD3.zmin)/2,['AR',num2str(ii),':',time])
		drawnow

		if (mod(ii,Skip)==0)
			%	mov(count) = getframe(gcf);
			DistanceMap = map.F;
			saveas(gcf, fullfile(Result_Folder,[num2str(ii),'AR','.png']))
			save(fullfile(Result_Folder,['DFV',num2str(ii),'AR','.mat']),'DistanceMap')
			%	count = count + 1;
		end
	end

	
	%saveas(gcf, fullfile(Result_Folder,[num2str(ii),'AR','.png']))

	

end

% write test end time

	Test_info = fopen(fullfile(Result_Folder,'test_end_info'), 'w');
	fprintf(Test_info, 'test end time : \t');
	fprintf(Test_info, [datestr(datetime('now'), 'yy/mm/dd HH:MM:SS'),'\n']);
	fclose(Test_info);

save(fullfile(Result_Folder,'movie.mat'),'mov')

	

