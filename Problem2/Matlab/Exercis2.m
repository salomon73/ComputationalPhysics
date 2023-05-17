%% Second exercise %% 

% Generate A and b 
    N=7;
    [A,b]=generate_ab(N); detA=det(A);
% Rank theorem
    rang=rank(A); 
    dimker=size(A,1)-rang;
    v1=null(A); 
    verif=norm(A*v1);
% Projecting on orthogonal complement of null(A)
    Proj=projection(N);
% verify non-singularity
    projA=Proj*A*Proj'; 
    detProjA=det(projA);
    y=projA\(Proj*b); 
    x=Proj'*y;
% Potential and current potted 
% through reshaping function
    [V,Ir,It]=reshaping(N,x);
% Kirschoff's law verification
    verif_kirschoff(Ir,It);
% Convergence study
% before running these lines, comment on the lines that create the figures 
% in the function reshaping 
    Req=V((N-1)/2,(N-1)/2)-V((N-1)/2+1,(N-1)/2+2);
    [Nmax,pi_num]=convergence(3,59);
    error=abs(pi-pi_num)/pi;
