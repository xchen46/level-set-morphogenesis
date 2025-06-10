function S = apply_Neumann(S,M,N,mode)

    if nargin < 4
        mode = 'symmetric';
    end

    [tot, ttot] = size(S);


   if strcmp(mode, 'one_sided')
       for boundary_LR = [1, M+1, tot - 2*M + 1, tot - M + 1]
           if boundary_LR == 1 || boundary_LR == tot - M + 1
               idx = boundary_LR : boundary_LR + M - 1;
           else
               idx = boundary_LR+1 : boundary_LR + M - 2;
           end

           S(sub2ind([tot, ttot], idx, idx)) = S(sub2ind([tot, ttot], idx, idx))-1;
           S(sub2ind([tot, ttot], idx, idx + (2 * (boundary_LR <= M+1) - 1) * M)) = S(sub2ind([tot, ttot], idx, idx + (2 * (boundary_LR <= M+1) - 1) * M))+1;
       end

       for boundary_UD = [M+1, M+2, 2*M-1, 2*M]
           idx = boundary_UD :M: tot-M;
           S(sub2ind([tot, ttot], idx, idx)) = S(sub2ind([tot, ttot], idx, idx))-1;
           S(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1))) = S(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1)))+1;
       end

   end



    if strcmp(mode, 'symmetric')
        for boundary_LR = [1, M+1, tot - 2*M + 1, tot - M + 1]
            if boundary_LR == 1 || boundary_LR == tot - M + 1
               idx = boundary_LR : boundary_LR + M - 1;
           else
               idx = boundary_LR : boundary_LR + M - 1;
               k = sub2ind([tot, ttot], idx, idx);
               g = sub2ind([tot, ttot], idx, idx - (2 * (boundary_LR <= M+1) - 1) * M);
               S(k) = S(k)-1;
               S(g) = S(g)+1;
            end

            k = sub2ind([tot, ttot], idx, idx);
            kk = sub2ind([tot, ttot], idx, idx + (2 * (boundary_LR <= M+1) - 1) * M);
            S(k) = S(k)-1;
            S(kk) = S(kk)+1;
            
        end


        for boundary_UD = [M+1, M+2, 2*M-1, 2*M]

            idx = boundary_UD-M: M : tot;

            if boundary_UD == M+1 || boundary_UD == 2*M
            else
                S(sub2ind([tot, ttot], idx, idx)) = S(sub2ind([tot, ttot], idx, idx))-1;
                S(sub2ind([tot, ttot], idx, idx - (2 * (boundary_UD < 2*M-1) - 1))) = S(sub2ind([tot, ttot], idx, idx - (2 * (boundary_UD < 2*M-1) - 1)))+1;
            end

            S(sub2ind([tot, ttot], idx, idx)) = S(sub2ind([tot, ttot], idx, idx))-1;
            S(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1))) = S(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1)))+1;

        end


        inner_LR = {
            2*M+1:M:tot-3*M+1,                      % edge 1 (e.g., top inner indices)
            3*M:M:tot-2*M            % edge 2 (e.g., bottom inner indices)
            };

        for edge = 1:2
            idx = inner_LR{edge};

            S(sub2ind([tot, ttot], idx, idx)) = S(sub2ind([tot, ttot], idx, idx))-2;
            S(sub2ind([tot, ttot], idx, idx + M)) = S(sub2ind([tot, ttot], idx, idx + M))+1;
            S(sub2ind([tot, ttot], idx, idx - M)) = S(sub2ind([tot, ttot], idx, idx - M))+1;
        end


        inner_UD = {
            3:M-2,                      % edge 1 (e.g., top inner indices)
            tot-N+3 : tot-2            % edge 2 (e.g., bottom inner indices)
            };

        for edge = 1:2
            idx = inner_UD{edge};

            S(sub2ind([tot, ttot], idx, idx)) = S(sub2ind([tot, ttot], idx, idx))-2;
            S(sub2ind([tot, ttot], idx, idx + 1)) = S(sub2ind([tot, ttot], idx, idx + 1))+1;
            S(sub2ind([tot, ttot], idx, idx - 1)) = S(sub2ind([tot, ttot], idx, idx - 1))+1;
        end
       


         % Collect all boundary rows (top, bottom, left, right)
       top_edge    = 1:M;
       bottom_edge = M*(N-1)+1 : M*N;
       left_edge   = 1:M:N*M;
       right_edge  = M:M:N*M;

       % Combine unique indices
       boundary_rows = unique([top_edge, bottom_edge, left_edge, right_edge]);
       S(boundary_rows,:) = S(boundary_rows,:)/2;

    end
  


end

