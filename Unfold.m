function [X] = Unfold( X, dim, i )
X = reshape(shiftdim(X,i-1), dim(i), []);%Ïò×óÒÆ¶¯i_1Î» reshape96*30*40