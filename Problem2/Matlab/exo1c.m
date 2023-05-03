%% LU decomposition of an 'ill' matrix %%

% 'ill' matrix
A=[1  2  3 ;
   2  4  9 ;
   4 -3  1 ];

% LU decomposition without pivoting
% result gives inf values
[L,U] = lu_decomposition_nopiv(A)

% Matlab and lu with pivoting
% correct answer
[L_matlab, U_matlab, P_matlab] = lu(A)
[L_pivot U_pivot P_pivot] = lu_decomposition(A)

