classdef ReconRadial2D
    %RECONRADIAL2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ser % array of projections
        nProj % number of projections
        nFid % number of data points per projection
        kr % normalized kr, from 0 to pi
        rho
        rho_nocorrection
        rho2
    end
    
    methods
        function obj = ReconRadial2D(ser)
            obj.ser = ser;
            [obj.nProj, obj.nFid] = size(ser);
            
            obj.kr = (pi / obj.nFid) * ( 0:(obj.nFid-1) );
            obj.kr = obj.kr(:);
            
        end
        
        function [ kx, ky ] = traj_radial_2D( obj, nProj, kr )
            %
            % Calculate 2D radial projections, normalized.
            % Np: Number of projections
            % kr: normalized kr, from 0 to pi
            %
            
            phi = ( 0:(nProj - 1) ) * 2.0 * pi / nProj;
            kx = kron(kr,cos(phi));
            ky = kron(kr,sin(phi));
            kx = kx(:);
            ky = ky(:);
            
        end
    
        function [ wk1, wk2 ] = density_radial_2d( obj, nProj, kr, flag )
            % Calculating the density weighting factor for the 2D radial MRI sequence
            %
            % Np: Number of projections
            % 
            % flag: Density for the center
            %   if flag = 0, pi*wdkr^2
            %   if flag = 1, 1/Np;
            %   if flag = 2, 
            
            wdkr = mean(diff(kr));
            wk1 = 2 * pi * kr .* wdkr / nProj;
            wk2 = wk1;
            
            switch flag
                case 0
                    wk1(1) = pi*(wdkr/2)^2;
                    wk2(1) = pi*(wdkr/2)^2 / nProj;
                case 1
                    wk1(1) = 1 / nProj; % not working?
                    wk2(1) = 1 / nProj;
                case 2
                    wk1(1) = wk1(1) / nProj;
                    wk2(1) = wk2(1) / nProj;
            end
            
            wk1 = repmat(wk1, [1, nProj]);
            wk1 = wk1(:);
            
            wk2 = repmat(wk2, [1, nProj]);
            wk2 = wk2(:);
            
        end
        
        function obj = recon2D(obj)
            tolNUFFT = 1e-6;
            
            fid_DC = obj.nFid - 20;
            ser2 = obj.ser - repmat( mean(obj.ser(:,fid_DC:end),2) , [1, obj.nFid]);
            
            first = 1;
            Nfid = 448;

            s = obj.ser(:,first:(Nfid+first-1))';
            s = s(:);
            
            nx = 256;
            kr = pi*(0:(Nfid-1))'/ (Nfid);
 
            [ kx, ky ] = obj.traj_radial_2D( obj.nProj, kr );
            
            figure
            subplot(1,2,1)
            scatter3(kx,ky,abs(s),1)
            subplot(1,2,2)
            scatter3(kx,ky,log(abs(s)),1)
            
            flag = 0;
            [ wk1, wk2 ] = obj.density_radial_2d( obj.nProj, kr, flag );
            wk1 = wk1/max(wk1(:));
            wk2 = wk2/max(wk2(:));
            
            tic 
            rho2 = nufft2d1(obj.nProj * Nfid, kx, ky, wk2 .* s, -1, tolNUFFT, nx, nx);
            toc
            
            wk3 = nufft2d1(obj.nProj * Nfid, kx, ky, wk1, -1, tolNUFFT, nx, nx);
            wk4 = nufft2d1(obj.nProj * Nfid, kx, ky, wk2, -1, tolNUFFT, nx, nx);
            
            figure
            imagesc(real(abs(rho2)))
            axis image
            
            figure
            mesh(abs(rho2))
           
            rho3 = fft2(fftshift(rho2));
            wk5 = fft2(fftshift(wk3));
            wk6 = fft2(fftshift(wk4));
            
            % Deconvolution of the density weights
            
            wmax = 0.001 * abs(wk5(1,1));
            toss = find( abs(wk5) <= wmax );
            keep = find( abs(wk5) > wmax );
            
            wk7(toss) = 0.0;
            wk7(keep) = 1./wk5(keep);
            wk7 = reshape(wk7, [nx, nx]);
            
            wmax = 0.001 * abs(wk6(1,1));
            toss = find( abs(wk6) <= wmax );
            keep = find( abs(wk6) > wmax );
            
            wk8(toss) = 0.0;
            wk8(keep) = 1./wk6(keep);
            wk8 = reshape(wk8, [nx, nx]);
            
            rho4a = fftshift(ifft2(rho3./wk5));
            rho4b = fftshift(ifft2(rho3./wk6));
            rho4 = fftshift(ifft2(rho3.*wk7));
            rho4c = fftshift(ifft2(rho3.*wk8));

            figure
            subplot(1,2,1)
            mesh(abs(rho4))
            subplot(1,2,2)
            mesh(abs(rho4c))
            
            figure
            subplot(1,2,1)
            imagesc(abs(rho4))
            axis image
            subplot(1,2,2)
            imagesc(abs(rho4c))
            axis image
            
            obj.rho_nocorrection = rho2;
            obj.rho = rho4;
            obj.rho2 = rho4c;

        end
    end

end

