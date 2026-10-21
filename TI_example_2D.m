for width = 30:10:40;


Noise_level=0.0;

for kappa = 0.02;

x0 = kappa*0;
x1 = kappa*width/2;
y0 = kappa*0.0;
y1 = kappa*0.0;
z0 = 0.0;
z1 = 0.0;

steps = 100;

file_nameA = "2DAL";

file_nameA = append(file_nameA,string(width));

file_nameA = append(file_nameA,"Kappa");

file_name_baseA = append(file_nameA,string(kappa))



[H_per,H,exp_x,exp_y,X,Y] = ChernSystem(width, Noise_level,1);
X1 = kappa*X;X2 =kappa*Y;X3 = H;
line_of_triple(X1,X2,X3,x0,x1,y0,y1,z0,z1,steps,file_name_baseA)



end
end









