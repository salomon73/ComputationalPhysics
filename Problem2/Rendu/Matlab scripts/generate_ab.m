%% Generator of A and b %% 
%This code is intended to automatically generate the A matrix of the
%problem. To do this, the Ohm and Kirschoff equations and the edge conditions
%have been written in the $N=5$ and $N=7$ cases in order to find a pattern and 
%write this code. 

function [A,b]=generate_ab(N)
    if mod(N,2)==0
        error('N has to be odd.')
    else
        Nu=3*N*(N-1); Nc=Nu;
        A=zeros(Nc,Nu); b=zeros(Nc,1);
    %Creating the matrix A 
        %Ohm's law
        shift_lign=0; shiftIt=0; shiftV=0;
        tempIt=zeros(N,N-2); tempV=zeros(N*(N-2),2);
        for i=1:(N-1)^2+(N-2)*N
            if i<=N
                for j=1:N-2
                    tempIt(i,j)=ind2sub(size(A),i+2*N*(N-1)+shiftIt);
                    tempV(1+shiftV,[1,2])=ind2sub(size(A),[i+shiftV i+1+shiftV]);
                    shiftIt=shiftIt+1; shiftV=shiftV+1;
                end
                It=reshape(transpose(tempIt),[],1);
                V=reshape(transpose(tempV),[],1);
            end
            if i<=(N-1)^2
                [c1,~]=ind2sub(size(A),i);
                [c2,~]=ind2sub(size(A),i+N-1);
                [c3,~]=ind2sub(size(A),i+N*(N-1));
                A(i,c1)=1; A(i,c2)=-1; A(i,c3)=-1; 
            else
                A(i,[V(i-1-N*(N-2)+shift_lign) V(i-N*(N-2)+shift_lign)])=[1 -1]; A(i,It(i-1-N*(N-2)))=-1;
                shift_lign=shift_lign+1;
            end
        end
        %Kirchhof's law
        shiftIt=0; shift_lines=0; 
        shiftIr=0; tempIr=zeros((N-1),(N-1));
        templign=zeros(N*(N-2),2); 
        %here we rearrange the vector It to have coherent eq.
        if N==3
           perm=flip(It(1:2),1);
           It(1:2,:)=perm;
        else 
           perm=It(1:N-2,:);
           It(1:N-2,:)=It((N-2)^2+1:(N-1)*(N-2),:);
           It((N-2)^2+1:(N-1)*(N-2),:)=perm;
           bound_r=N-2;
           for i=1:N-3
               perm=It((i-1)*bound_r+1:i*bound_r,:);
               It((i-1)*bound_r+1:i*bound_r,:)=It(i*bound_r+1:(i+1)*bound_r,:);
               It(i*bound_r+1:(i+1)*bound_r,:)=perm;
           end
        end
        %implementation of Kirschoff's eq.
        for i=(N-1)^2+(N-2)*N+1:(N-1)^2+(N-2)*N+N*(N-1)    
            if i<=(N-1)^2+(N-2)*N+N
               for j=1:N-1
                   tempIr(i-((N-1)^2+(N-2)*N),j)=ind2sub(size(A),(N-1)*N+1+shiftIr);
                   shiftIr=shiftIr+1;
               end
               Ir=reshape(transpose(tempIr([1 (N-1)],:)),[],1);
               for j=1:N-2
                   templign(1+shiftIt,[1,2])=ind2sub(size(A),[i+shiftIt i+1+shiftIt]);
                   shiftIt=shiftIt+1;
               end
               lines=reshape(transpose(templign),[],1);
            end
            if i<=(N-1)^2+(N-2)*N+(N-1)*(N-2) 
                [c1,~]=ind2sub(size(A),N*(N-1)+i-((N-1)^2+(N-2)*N));
                [c2,~]=ind2sub(size(A),N*(N-1)+i-((N-1)^2+(N-2)*N)+N-1);
                A(i,[c1 c2])=[-1 1]; 
            else
                if i<=(N-1)^2+(N-2)*N+(N-1)*(N-2)+(N-1)
                   A(i,Ir(i-((N-1)^2+(N-2)*N+(N-1)*(N-2))))=1;
                else
                   A(i,Ir(i-((N-1)^2+(N-2)*N+(N-1)*(N-2))))=-1;
                end               
            end
            if i<=(N-1)^2+(N-2)*N+length(lines)/2
            A([lines(i-((N-1)^2+(N-2)*N)+shift_lines) lines(i-((N-1)^2+(N-2)*N)+1+shift_lines)],It(i-((N-1)^2+(N-2)*N)))=[1 -1];
            shift_lines=shift_lines+1;           
            end
        end
        %BC
        tempIr=(N-1)*N+(N-1)^2+1:1:2*(N-1)*N;
        tempIt=2*N*(N-1)+N-1:N-1:Nu;
        Ir=reshape(transpose(tempIr),[],1);
        It=reshape(transpose(tempIt),[],1);
        for i=(N-1)^2+(N-2)*N+N*(N-1)+1:Nc
            if i<=(N-1)^2+(N-2)*N+N*(N-1)+N-1
                A(i,Ir(i-((N-1)^2+(N-2)*N+N*(N-1))))=1;
            else 
                A(i,It(i-((N-1)^2+(N-2)*N+N*(N-1)+N-1)))=1;
            end
        end
    %rearranging the matrix A to have coherent eq.
        Ohm=A(1:(N-1)^2+(N-2)*N,:);
        Kirschoff1=A((N-1)^2+(N-2)*N+1:(N-1)^2+(N-2)*N+(N-2)*(N-1),:);
        Kirschoff2=A((N-1)^2+(N-2)*N+(N-2)*(N-1)+1:(N-1)^2+(N-2)*N+(N-2)*(N-1)+N-1,:);
        Kirschoff3=A((N-1)^2+(N-2)*N+(N-2)*(N-1)+N:(N-1)^2+(N-2)*N+N*(N-1),:);
        A(1:N-1,:)=Kirschoff2;
        A(N:N+(N-2)*(N-1)-1,:)=Kirschoff1;
        A(N+(N-2)*(N-1):N*(N-1),:)=Kirschoff3;
        A(N*(N-1)+1:N*(N-1)+(N-1)^2+(N-2)*N,:)=Ohm;
    %Creating vector b
        bin=sub2ind([N-1 N],(N-1)/2,(N-1)/2);
        bout=sub2ind([N-1 N],(N-1)/2+1,(N-1)/2+2);
        b([bin bout])=[1 -1];
    end
end