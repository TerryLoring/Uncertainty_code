%
% 

function line_of_triple(X1,X2,X3,x0,x1,y0,y1,z0,z1,Nsteps,file_name_base)



X1_bracket_X2 =   i*(X1*X2 - X2*X1);
X1_bracket_X3 =   i*(X1*X3 - X3*X1);
X2_bracket_X3 =   i*(X2*X3 - X3*X2);
X1_bracket_X2 = (1/2)*(X1_bracket_X2'+X1_bracket_X2);
X1_bracket_X3 = (1/2)*(X1_bracket_X3'+X1_bracket_X3);
X2_bracket_X3 = (1/2)*(X2_bracket_X3'+X2_bracket_X3);



%
file_name1 = strcat(file_name_base,'_1.pdf');
file_name2 = strcat(file_name_base,'_2.pdf');


% The function eigs() throws lots of errors we want
% to ignore, so I make them be warngings instead.
warning('error','MATLAB:eigs:NoEigsConverged');
warning('error','MATLAB:nearlySingularMatrix');
warning('error','MATLAB:eigs:SigmaNearExactEig');
warning('error','MATLAB:eigs:IllConditionedA');

n = size(X1);
n = n(2);

I2 = sparse([1 0;0 1]);
sigma1 = sparse([0 1;1 0]);
sigma2 = sparse([0 -i;i 0]);
sigma3 = sparse([1 0;0 -1]);

g1 = sigma1;
g2 = sigma2;
g3 = sigma3;

X1_squared = X1*X1;
X2_squared = X2*X2;
X3_squared = X3*X3;
Q0 = X1_squared + X2_squared + X3_squared;



I_n = speye(n);

	
function Q = quad_composite(x1,x2,x3)
	Q = Q0 + (-2*x1)*X1 + (-2*x2)*X2 + (-2*x3)*X3 + (x1^2+x2^2+x3^2)*I_n;
	Q = (1/2)*(Q' + Q);
end

function LL = localizer_squared(x1,x2,x3)
	LL = kron(I2,quad_composite(x1,x2,x3));
	LL = LL + kron(-i*g1*g2,X1_bracket_X2);
	LL = LL + kron(-i*g1*g3,X1_bracket_X3);
	LL = LL + kron(-i*g2*g3,X2_bracket_X3);
	LL = (1/2)*(LL' + LL);
end


function u = uncert(A,v)
	u = abs(v'*(A^2)*v - (v'*A*v)^2);
end

function u = expectation(A,v)
	u = real(v'*(A)*v);
end



%Reserve memory for data arrays
t_values = zeros(size(1,Nsteps));
Cliff_gap = zeros(size(1,Nsteps));
quad_gap = zeros(size(1,Nsteps));
sum_uncertainty = zeros(size(1,Nsteps));
RMS_expectation = zeros(size(1,Nsteps));

for j = 1:Nsteps
	t = (j-1)/(Nsteps-1);
	t_values(j)=t;
	x = (1-t)*x0+t*x1;
	y = (1-t)*y0+t*y1;
	z = (1-t)*z0+t*z1;
	
	LL = localizer_squared(x,y,z);
	try
		opts.tol = 10^-14;
		opts.p = 60;
		opts.maxit = 600;
			e_val = eigs(LL,1,'sm',opts);
			gap_Cliff = sqrt(abs(e_val));
	catch
		gap_Cliff = 0;
	end	
	Cliff_gap(j) = gap_Cliff;
	
	
	Q = quad_composite(x,y,z);
	try
		opts.tol = 10^-14;
		opts.p = 60;
		opts.maxit = 500;
			[e_vec,e_val] = eigs(Q,1,'sm',opts);
			gap_quad = sqrt(abs(e_val));
	catch
		gap_quad = 0;
		e_vec = (1/sqrt(n))*ones(n,1);
	end	
	quad_gap(j) = gap_quad;
	
	
	sum_uncert = uncert(X1,e_vec) ...
	           + uncert(X2,e_vec) ...
	           + uncert(X3,e_vec);
	sum_uncertainty(j) = sum_uncert;
	
	RMS_exp = (expectation(X1_bracket_X2,e_vec))^2 ...
	        + (expectation(X1_bracket_X3,e_vec))^2 ...
	        + (expectation(X2_bracket_X3,e_vec))^2;
	RMS_exp = sqrt(RMS_exp);
	RMS_expectation(j) = RMS_exp;
	
end


figure('Visible', 'off')
hold on
plot(t_values,quad_gap, 'k');
plot(t_values,Cliff_gap,'k--');
plot(t_values,sqrt(Cliff_gap.^2+RMS_expectation),'k:');
currentLimits = ylim;
ylim([0 currentLimits(2)])
drawnow
exportgraphics(gcf,file_name1)
pause(0.1)


figure('Visible', 'off')
hold on
plot(t_values,quad_gap,'k')
plot(t_values,sqrt(sum_uncertainty),'k--')
currentLimits = ylim;
ylim([0 currentLimits(2)])
drawnow
exportgraphics(gcf,file_name2)

end
