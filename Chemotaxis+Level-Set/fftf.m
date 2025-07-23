close all
clc
%%
membrane_size = 50;
grid_size = 256;
time_end = size(U,2);
dx = membrane_size/(grid_size - 1);  % Grid step

predicted_length = 2*pi/1.4651;
area_of_interests = ceil(predicted_length*4/dx);

r = randi([area_of_interests+3 grid_size-3-area_of_interests],1,2);

n = area_of_interests + 1;
L = (n-1)*dx;
w = tukeywin(n, 0.25);         % or hann(n)
padF = 4;
pad = padF * n;
W = w * w';


z = reshape(U(:,time_end),grid_size,grid_size);
z = z(r(1):r(1)+area_of_interests,r(2):r(2)+area_of_interests);
z0 = (z - mean(z(:))) .* W;


Uhat = fftshift(fft2(z0,pad,pad));
spec = abs(Uhat).^2;

%%
dk       = 2*pi/(L*padF);               % one FFT step
k  = (-pad/2 : pad/2-1) * dk;  % same dk as above
[KX,KY] = meshgrid(k,k);
K  = sqrt(KX.^2 + KY.^2);

%%
Kedges   = 0 : dk : max(K(:))/5;       % MANY narrow annuli
Kbins = 0.5*(Kedges(1:end-1)+Kedges(2:end));
spec_r = zeros(size(Kbins));

for i = 1:length(Kedges)-1
    mask = (K >= Kedges(i)) & (K < Kedges(i+1));
    spec_r(i) = mean(spec(mask), 'all');
end
%%
[~, idx] = max(spec_r);   % ignore DC (k=0)
k_star      = Kbins(idx);
lambda_char = 2*pi ./ Kbins;
fprintf('Δk = %.4f   k* = %.4f   λ* = %.4f\n', dk, k_star, 2*pi/k_star);
%%

% 2D spectrum
imagesc(k, k, log1p(spec));
axis equal tight; colorbar;
title('Log Power Spectrum'); xlabel('k_x'); ylabel('k_y');

% 1D radially averaged
figure;
plot(Kbins, spec_r);
xlabel('wavenumber k'); ylabel('power');
title('Radial average of 2D spectrum');


%%

% for i = 1:1000
%     if mod(i,10) == 0
%         clf
% 
%              % power spectrum
% 
% 
%         % imagesc(z(:,end-1:end)), colorbar;
%         imagesc(spec); axis equal tight; colorbar
% 
%         pause(0.1)
%         title('Concentration of U');
% 
%         % frame = getframe(gcf);
%         % writeVideo(vid, frame);
% 
%     end
% end 