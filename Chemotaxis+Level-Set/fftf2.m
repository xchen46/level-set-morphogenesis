%%
clc
close all
load('comparsion1_c2.mat');
%%

membrane_size = 50;
grid_size = 250;
time_end = size(U,2);
dx = membrane_size/(grid_size - 1);  % Grid step

% k_theory = 1.55891 ;
k_theory = 1.4651 ;
predicted_length = 2*pi/k_theory;
area_of_interests = ceil(predicted_length*4/dx);
z = reshape(U(:,time_end),grid_size,grid_size);
%%
n       = area_of_interests + 1;      % ← keep your existing size
L       = (n-1)*dx;
stride  = floor(n/2);                 % 50 % overlap (= smoother map)
margin  = 4;                          % avoid 3-px boundary artefacts

w = tukeywin(n, 0.25);         % or hann(n)
padF = 4;
pad = padF * n;

%%

ix  = margin          : stride : (grid_size - margin - n);
iy  = ix;                                 % square scan
kMap = NaN(numel(ix), numel(iy));         % local k* field

%%

% --- fixed stuff you already compute once ---
dk      = 2*pi/(L*padF);
[KX,KY] = meshgrid((-pad/2:pad/2-1)*dk);  K = sqrt(KX.^2+KY.^2);
Kedges  = 0 : dk : 5;             % full radial range
Kbins   = 0.5*(Kedges(1:end-1)+Kedges(2:end));
W2D     = w * w.';                        % reuse Tukey window


for a = 1:numel(ix)
  for b = 1:numel(iy)
      % ---- crop one window -------------------------------------------
      x0 = ix(a);  y0 = iy(b);
      patch = z( x0:x0+n-1 , y0:y0+n-1 );
      patch = (patch-mean(patch(:))).*W2D;    % window & demean

      % ---- FFT & power spectrum --------------------------------------
      spec = abs( fftshift( fft2(patch, pad, pad) ) ).^2;


      [~, idxMax] = max(spec(:));                  % brightest pixel in 2-D
      [row,col]   = ind2sub(size(spec), idxMax);
      kMap(a,b)   = hypot(KX(row,col), KY(row,col));   % physical k*
      % % % ---- radial average & local k* ---------------------------------
      % spec_r = zeros(size(Kbins));
      % for j = 1:numel(Kbins)
      %     mask = (K>=Kedges(j)) & (K<Kedges(j+1));
      %     spec_r(j) = mean(spec(mask),'all');
      % end
      % 
      % [~, iPeak] = max(spec_r);
      % kMap(a,b)  = Kbins(iPeak); 
  end
end
%%

figure, imagesc(iy, ix, kMap), axis equal tight
colorbar, title('local k* map (rad / unit length)')

kAll  = kMap(~isnan(kMap));        % vector of all local peaks
k_star = mean(kAll);
err_pct  = 100*(k_star - k_theory)/k_theory;
fprintf('Mean k* = %.4f   SD = %.4f   Err = %.4f\n', ...
        mean(kAll), std(kAll), err_pct);
