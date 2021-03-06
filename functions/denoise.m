function [Z,s] = denoise(X,n)
%PCA denoising of 2D relaxation data
%Input is mxn matrix, with m spectral points and n relaxation points
%n is the number of principal components to retain after denoising

%Returns Z denoised matrix and s singular values vector

%Removed dependence on Machine Leearning Toolbox for PCA
%[V, U] = pca(X,'Algorithm','svd');   %U will be nechoes x nechoes -> V is PC coeff, U is PC

[~,b] = size(X);
for j = 1:b %'center' the matrix to have mean 0 on the columns
    X(:,j) = X(:,j) - mean(X(:,j)) ;
end

[U, s, V] = svd(X); %singular values, U SIG V 
B = s*V';
U = U(:,1:n);
B = B(1:n,:);
Z=U*B;
s=diag(s);
end

