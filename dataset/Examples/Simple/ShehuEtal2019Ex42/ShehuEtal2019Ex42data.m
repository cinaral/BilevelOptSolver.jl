function data = ShehuEtal2019Ex42data(row,ny,mu)
 rng('default'); rng(1);
P      = randn(ny,ny);
data.Q = P*P'/ny; 
%data.Q = eye(ny);

k      = floor(0.05*ny);
y      = zeros(ny,1);
I0     = randperm(ny); 
I      = I0(1:k);
y(I)   = rand(k,1)+5; % generate a sparse vector y

data.A = randn(row,ny);
data.b = data.A(:,I)*y(I)+0.01*randn(row,1);
data.y = y; 
data.mu= mu;
data.At= data.A';
data.AA= data.At*data.A;
end

