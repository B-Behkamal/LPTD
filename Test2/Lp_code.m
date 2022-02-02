function [permvec,FVAL] = Lp_code(D,min)

    [N1,N2] = size(D);
    if N1 <= N2 
        [permvec,res,FVAL] = injmatch(D,min);
    else
        [permvec,res,FVAL] = injmatch(D',min);
        permvecold = permvec;
        permvec = zeros(N1,1);
        for i=1:N2
            permvec(permvecold(i))= i;
        end
    end
end



function [permvec,res,FVAL]=injmatch(D, min)  
    [N1,N2] = size(D);
    assert(N1<=N2)
    
    if (min == 0)
        vD = D(:);
    else
        vD= -D(:);
    end


    lb = zeros(N1*N2,1);
    ub = ones(N1*N2,1);

    Aeq = zeros(N1,N1*N2);
    beq = ones(N1,1);    
    for i=1:N1
        vecidx = i:N1:N1*N2;
        Aeq(i,vecidx) = 1;
    end
    
    A = zeros(N2,N1*N2);
    b = ones(N2,1); 
    offset = 0;
    for i=1:N2        
        for j=1:N1
            A(i,j+offset) = 1;
        end
        offset = offset + N1;
    end
    
    [res,FVAL]=linprog(vD,A,b,Aeq,beq,lb,ub);
    permvec = generate_perm(res,N1,N2);

end



function permvec=generate_perm(res,N1,N2)
    
    mres = reshape(res,N1,N2);    % matrix result (row=stick and column=helix)
    permvec = zeros(N1,1);
    ctr = 1;
    for i=1:N1
        for j=1:N2
            if(mres(i,j) > 0.5)
                permvec(ctr) = j;
                ctr = ctr + 1;
                break
            end
        end
    end 
    if ctr ~= N1+1
        warning('not a permutation');
    end
end