% builds the Hamiltonians and position observables for both periodic
% and zero boundary conditions.


% The Hamiltonian is the one Hastings and I used, just the spin up only
% part with  Delta=0,  in our EPL paper.
%
% --Terry A. Loring, February 11,2015.


function [H_per,H,exp_x,exp_y,X,Y] = ChernSystem(width, Noise_level,seed_for_rand)

% Should be A=1; B=-1; C=0; D=0; M =-1 to recreate Figure 1 in EPL paper.
% was A=1; B=-1; C=0; D=0; M =-1; (in annals of physics paper and some parts of EPL paper)

A=1; B=-1; C=0; D=0; M =-2;

% A=1 here corresponds to A=2 in standard picture....

bandE=0;
bandH=1;

height = width;

dim = 2*width*height;


% To convert from a position, band and spin to matrix row or column index
function k = ind(x,y,band)
	w1 = width;
	w2 = w1*height;
	w3 = w2*2;
	k = 1 + x + y*w1 + band*w2;
end


%Build the Lattice, diagonal position observables
exp_x = sparse(dim,dim);
exp_y = sparse(dim,dim);
X = sparse(dim,dim);
Y = sparse(dim,dim);
for Band = 0:1
	for ind_x = 0:(width-1)
		for ind_y = 0:(height-1)
			k = ind(ind_x,ind_y,Band);
			exp_x(k,k) = exp(2*pi*i*ind_x/width);
			exp_y(k,k) = exp(2*pi*i*ind_y/height);
			X(k,k) = ind_x-width/2 + 1/2;
			Y(k,k) = ind_y-height/2 + 1/2;
		end
	end
end

rng(seed_for_rand);


%Build the Hamiltonian
%H = sparse(dim,dim);
Iall = zeros(1,0);
Jall = zeros(1,0);
Sall = zeros(1,0);
for ind_x = 0:(width-1)
	for ind_y = 0:(height-1)
	
		%Put in terms for staying at a site
		epsilon_k = C - 4*D;
		noise = Noise_level*(rand-0.5);
		%noise = normrnd(0,Noise_level);
		stationary = epsilon_k + noise;
		Mdiag = M-4*B;
		for Band = 0:1
			j = ind(ind_x,ind_y,Band);
			%H(j,j) = stationary + Mdiag*(2*Band-1);
			Iall(end+1) = j; Jall(end+1) = j;
			Sall(end+1) = stationary + Mdiag*(2*Band-1);
		end

		
		%Put in terms for North-South Neighbors
		north_x = mod(ind_x+1,width);
		for Band=0:1
				j = ind(ind_x,ind_y,Band);
				k = ind(north_x,ind_y,Band);
				%H(j,k) = D+B*(2*Band-1);
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = D+B*(2*Band-1);
				j = ind(ind_x,ind_y,Band);
				k = ind(north_x,ind_y,1-Band);
				%H(j,k) = A*i;
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = A*i;
				j = ind(north_x,ind_y,Band);
				k = ind(ind_x,ind_y,Band);
				%H(j,k) = D+B*(2*Band-1);
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = D+B*(2*Band-1);
				j = ind(north_x,ind_y,Band);
				k = ind(ind_x,ind_y,1-Band);
				%H(j,k) = -A*i;
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = -A*i;
		end
		
		%Put in terms for East-West Neighbors
		east_y = mod(ind_y+1,height);
		for Band=0:1
				j = ind(ind_x,ind_y,Band);
				k = ind(ind_x,east_y,Band);
				%H(j,k) = D+B*(2*Band-1);
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = D+B*(2*Band-1);
				j = ind(ind_x,ind_y,Band);
				k = ind(ind_x,east_y,1-Band);
				%H(j,k) = A*(2*Band-1);
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = A*(2*Band-1);
				j = ind(ind_x,east_y,Band);
				k = ind(ind_x,ind_y,Band);
				%H(j,k) =  D+B*(2*Band-1);
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = D+B*(2*Band-1);
				j = ind(ind_x,east_y,Band);
				k = ind(ind_x,ind_y,1-Band);
				%H(j,k) = -A*(2*Band-1);
				Iall(end+1) = j; Jall(end+1) = k;
				Sall(end+1) = -A*(2*Band-1);
		end
		
	end	
end
H = sparse(Iall,Jall,Sall,dim,dim);

H_per = H;

% Zero as boundary conditions
ind_x = 0;
north_x = width-1; 
for ind_y = 0:(height-1)
	for Band=0:1
		for altB=0:1
			j = ind(ind_x,ind_y,Band);
			k = ind(north_x,ind_y,altB);	
			H(j,k)=0;
			H(k,j)=0;
		end
	end
end
ind_y = 0;
east_y = height-1; 
for ind_x = 0:(width-1)
	for Band=0:1
		for altB=0:1
			j = ind(ind_x,ind_y,Band);
			k = ind(ind_x,east_y,altB);	
			H(j,k)=0;
			H(k,j)=0;
		end
	end
end


end
