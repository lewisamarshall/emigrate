function [A,B]=PrepareDerivMatrices(x,Order,Flag)

% Implicit, 6th order compact calculation of first derivative
% A*f'=B*f

h=x(2)-x(1);
% Assuming uniform grid spacing
N=size(x,1);
A=sparse(N,N);    B=sparse(N,N);

switch (Flag)
    case ('Compact')
        if Order==1
            % Forward derivative for first(/last) grid point - 5th order
            A_First=[1 4];
            B_First=(1/h)*[-37/12    2/3    3   -2/3    1/12];

            % Skewed derivative for second/(N-1) grid point - 5th order
            A_Second=[1/6 1 1/2];
            B_Second=(1/h)*[-10/18   -1/2    1    1/18];


            % Central derivative for internal grid points
            alphaI=1/3;
            A_Internal=[alphaI 1 alphaI];
            B_Internal=[-1/3*(4*alphaI-1)*(1/(4*h))  -2/3*(alphaI+2)*(1/(2*h))  0  2/3*(alphaI+2)*(1/(2*h))  1/3*(4*alphaI-1)*(1/(4*h))];


            A(1,1:length(A_First))      = A_First;        B(1,1:length(B_First)) = B_First;
            A(2,1:length(A_Second))     = A_Second;       B(2,1:length(B_Second))= B_Second;

            A(end,end:-1:end-length(A_First)+1)    = A_First;       B(end,end:-1:end-length(B_First)+1) = - B_First;
            A(end-1,end:-1:end-length(A_Second)+1) = A_Second;      B(end-1,end:-1:end-length(B_Second)+1) = - B_Second;


            for ij=3:N-2  % run on all rows
                A(ij,ij-1:ij+1)=A_Internal;
                B(ij,ij-2:ij+2)=B_Internal;
            end

        elseif Order==2
            % Forward derivative for first(/last) grid point - 5th order
            A_First=[1 137/13];
            B_First=(1/h^2)*[1955/156 -4057/156 1117/78 -55/78 -29/156 7/156];

            % Skewed derivative for second/(N-1) grid point - 5th order
            A_Second=[1/10 1 -7/20];
            B_Second=(1/h^2)*[99/80   -3    93/40   -3/5    3/80];


            % Central derivative for internal grid points
            alphaI=2/11;
            aI=4/3*(1-alphaI)*(1/h^2);
            bI=1/3*(-1+10*alphaI)*(1/(4*h^2));
            A_Internal=[alphaI 1 alphaI];
            B_Internal=[bI   aI  -2*(aI+bI)  aI  bI];


            A(1,1:length(A_First))      = A_First;        B(1,1:length(B_First)) = B_First;
            A(2,1:length(A_Second))     = A_Second;       B(2,1:length(B_Second))= B_Second;

            A(end,end:-1:end-length(A_First)+1)    = A_First;       B(end,end:-1:end-length(B_First)+1) = B_First;
            A(end-1,end:-1:end-length(A_Second)+1) = A_Second;      B(end-1,end:-1:end-length(B_Second)+1) = B_Second;


            for ij=3:N-2  % run on all rows
                A(ij,ij-1:ij+1)=A_Internal;
                B(ij,ij-2:ij+2)=B_Internal;
            end

        end % if

    case {'2nd','Upwind','SLIP'}
        if Order==1
            A=eye(N);
            B=diag(ones(1,N-1)/(2*h),1)+diag(-ones(1,N-1)/(2*h),-1);
            B(1,1:2)=[-1/h,1/h];
            B(end,end-1:end)=[-1/h,1/h];
        elseif Order==2
            A=eye(N);
            B=diag(ones(1,N-1)/h^2,1)+diag(ones(1,N-1)/h^2,-1)-2*diag(ones(1,N)/h^2,0);
            B(1,1:3)=[1,-2,1]/h^2;
            B(end,end-2:end)=[1,-2,1]/h^2;
        end
end % switch
