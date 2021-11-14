function [Z,s] = denoise(spec,n)
%PCA denoising of 2D relaxation data
%Input is mxn matrix, with m spectral points and n relaxation points
%n is the number of principal components to retain after denoising

%Returns Z denoised matrix
%Returns s singular values

X  = spec;
[U1, s, V1] = svd(X); %singular values, U SIG V 
[V, U] = pca(X,'Algorithm','svd');   %U will be nechoes x nechoes -> V is PC coeff, U is PC
U = U(:,1:n);
V = V(:,1:n);
Z = U*V';
end

