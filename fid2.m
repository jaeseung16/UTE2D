gyro = 2 * pi * 42.6; % 42.6 MHz/T
R = 0.0025; % sample radius = 2.5 mm

N = 512; % data size

dw = 5; % dwell time in us

t = (0:N-1) * dw;

G = 0.03 * 0.52; % Graident in T/m

gammaG = gyro * G;

%%

[X Y] = meshgrid(-2:2, -2:2);

location = [0.015 * sqrt(X(:).^2 + Y(:).^2), atan2(Y(:), X(:))];

%%

omega = repmat( 2 * pi * [-0.0014; -0.0007; 0; 0.0007; 0.0014], [5, 1] );
radius = R * ones( length(X(:)), 1 );
density = ones( length(X(:)), 1 );
T2 = repmat( [1000, 100, 10, 1, 0.1] * 1000, [5, 1]);
T2 = T2(:);

nProj=512;

%%
[ fid ] = sumOverCircles( location, density, radius, T2, omega, gammaG, t, nProj );
UTE2D = ReconRadial2D( fid' );
UTE2D = UTE2D.recon2D

%%
n = 256;
dk = gyro * G * dw * 448 / n;
dx = pi / dk / n;
dy = dx;

original = zeros(n,n);

coordX = X(:) * 0.015;
coordY = Y(:) * 0.015;

for k = 1:n
    x = (k - n/2) * dx - dx;
    
    for l = 1:n
        y = (l - n/2) * dy - dy;
        
        for m = 1:length(coordX)
            if (coordX(m) - x)^2 + (coordY(m) - y)^2 <= R^2
                original(k,l) = 1;
            end
        end
       
    end
end

            
%%
figure
subplot(2,2,1)
imagesc(original)
axis image
subplot(2,2,2)
imagesc(abs(UTE2D.rho))
axis image
subplot(2,2,3)
mesh(original)
subplot(2,2,4)
mesh(abs(UTE2D.rho))