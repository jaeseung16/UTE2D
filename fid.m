gyro = 2 * pi * 42.6; % 42.6 MHz/T
R = 0.0025; % sample radius = 2.5 mm

N = 512; % data size

dw = 5; % dwell time in us

t = (0:N-1) * dw;

G = 0.1 * 0.52; % Graident in T/m

%%
kR = gyro * G * R * t;

% J1 = besselj(1, kR);
% 
% fid1 = 2.0 * pi * R^2 * J1 ./ kR;
% 
% fid1(1) = pi * R^2; % At t = 0, J1 / kR = 0.5

J0 = besselj(0, kR);
J2 = besselj(2, kR);
fid1 = 2.0 * pi * R^2 * ( J0 + J2 ) / 2.0; % 2 * J1(x) = x * ( J0(x) + J2(x) )

ser1 = repmat( fid1, [ N, 1 ] );

%%
% phase shift due to the location
r0 = 0.0075;
theta0 = (pi / 180.0) * 90.0 + 2.0 * pi * (0:N-1) / N;
kr0 = gyro * G * r0 * cos( theta0 )' * t;
phase = exp(1i * kr0);

ser2 = phase .* ser1;

%%

T2 = 100; % T2 in ms
decay = repmat( exp( - t / (T2 * 1000) ), [ N, 1 ]);

ser = (ser1 + ser2) .* decay .* repmat( exp( 1i * 2 * pi * 0.00 * t), [N 1] ) ;

plot( t / 1000, real(ser(1,:)), t / 1000, imag(ser(1,:)) )

%%
UTE2D = ReconRadial2D( ser )

%%

UTE2D = UTE2D.recon2D

%%
n = 256;
dk = gyro * G * dw * 448 / n;
dx = pi / dk / n;
dy = dx;

original = zeros(n,n);

for k = 1:n
    x = (k - n/2) * dx - dx;
    
    for l = 1:n
        y = (l - n/2) * dy - dy;
        
        if sqrt(x^2 + y^2) <= R
            original(k,l) = 1;
        end
        
        if sqrt( (x-r0*cos((pi / 180.0) * 90.0))^2 + (y-r0*sin((pi / 180.0) * 90.0))^2 ) <= R
            original(k,l) = 1;
        end
       
    end
end

figure
subplot(2,2,1)
imagesc(original)
axis image
subplot(2,2,2)
imagesc(abs(UTE2D.rho))
axis image
subplot(2,2,3)
imagesc(abs(UTE2D.rho2))
axis image
subplot(2,2,4)
imagesc(original - abs(UTE2D.rho) / max( abs( UTE2D.rho(:) ) ) )
axis image

figure 
subplot(1,3,1)
mesh(original)
subplot(1,3,2)
mesh(abs(UTE2D.rho) / max( abs( UTE2D.rho(:) ) ) )
subplot(1,3,3)
mesh(original - abs(UTE2D.rho) / max( abs( UTE2D.rho(:) ) ) )

%%
reconstructed = abs(UTE2D.rho);

psf = fftshift( ifft2( fft2(abs(UTE2D.rho)) ./ fft2(original) ) );

mesh(abs(psf))

%%

mesh( abs(conv2( original, psf, 'same' ) ) )
