% This script finds the solution for the exercise 1.1 linear system
% using the matlab method for matrix quotient: A/B

A=[ 2 1 -1  5 ;
    1 2  3 -1 ;
    1 0  1  6 ;
    1 3 -1  5 ];

B = [13; 37; 30; 19];

X=A\B;

% Solution : [2; 4; 10; 3;]
