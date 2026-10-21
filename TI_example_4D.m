for width = [5,8];

file_name = "4DAedgeL";
file_name = append(file_name,string(width));
file_name = append(file_name,"Kappa");

steps = 100;

%% init:
nUC = width;
t = 1;
t1 = 0.8*exp(1i*pi*0.1);
t2 = 0.0;
m = 0.5;
a = 1;

[HH,XX,YY,ZZ,WW] = ai_4d_hamil_pos(nUC,nUC,nUC,nUC,t,t1,t2,m,a);

kappa = 0.10;
	x0 = kappa*0;
	x1 = kappa*width/2;
	y0 = kappa*0;
	y1 = kappa*0;
	z0 = kappa*0;
	z1 = kappa*0;
	w0 = kappa*0;
	w1 = kappa*0;
	e0 = 0.0;
	e1 = 0.0;
	
	file_name_base = append(file_name,string(kappa))
	X1 = kappa*XX; X2 =kappa*YY; X3 =kappa*ZZ; X4 = kappa*WW; X5=HH;
	line_of_quint(X1,X2,X3,X4,X5,x0,x1,y0,y1,z0,z1,w0,w1,e0,e1,steps,file_name_base)


end