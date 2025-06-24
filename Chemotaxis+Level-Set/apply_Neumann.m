function [S1, S2] = apply_Neumann(S,M,N,mode)

    if nargin < 4
        mode = 'symmetric';
    end

    [tot, ttot] = size(S);
    S1 = S;

   if strcmp(mode, 'one_sided')
       S2 = S;
       % for boundary_LR = [1, M+1, tot - 2*M + 1, tot - M + 1]
       for boundary_LR = [1, tot - M + 1]

           if boundary_LR == 1 || boundary_LR == tot - M + 1
               idx = boundary_LR : boundary_LR + M - 1;
           else
               idx = boundary_LR+1 : boundary_LR + M - 2;
           end

           S1(sub2ind([tot, ttot], idx, idx)) = S1(sub2ind([tot, ttot], idx, idx))-1;
           S1(sub2ind([tot, ttot], idx, idx + (2 * (boundary_LR <= M+1) - 1) * M)) = S1(sub2ind([tot, ttot], idx, idx + (2 * (boundary_LR <= M+1) - 1) * M))+1;
       end

       % for boundary_UD = [M+1, M+2, 2*M-1, 2*M]
       for boundary_UD = [M+1, 2*M]
           idx = boundary_UD :M: tot-M;
           S2(sub2ind([tot, ttot], idx, idx)) = S2(sub2ind([tot, ttot], idx, idx))-1;
           S2(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1))) = S2(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1)))+1;
       end

   end



    if strcmp(mode, 'symmetric')
        for boundary_LR = [1, M+1, tot - 2*M + 1, tot - M + 1]
        % for boundary_LR = [1, tot - 2*M + 1]

            if boundary_LR == 1 || boundary_LR == tot - M + 1
               idx = boundary_LR : boundary_LR + M - 1;
           else
               idx = boundary_LR : boundary_LR + M - 1;
               k = sub2ind([tot, ttot], idx, idx);
               g = sub2ind([tot, ttot], idx, idx - (2 * (boundary_LR <= M+1) - 1) * M);
               S1(k) = S1(k)-1;
               S1(g) = S1(g)+1;
            end

            k = sub2ind([tot, ttot], idx, idx);
            kk = sub2ind([tot, ttot], idx, idx + (2 * (boundary_LR <= M+1) - 1) * M);
            S1(k) = S1(k)-1;
            S1(kk) = S1(kk)+1;
            
        end


       for boundary_UD = [M+1, M+2, 2*M-1, 2*M]
       % for boundary_UD = [M+1, 2*M]

            idx = boundary_UD-M: M : tot;

            if boundary_UD == M+1 || boundary_UD == 2*M
            else
                S1(sub2ind([tot, ttot], idx, idx)) = S1(sub2ind([tot, ttot], idx, idx))-1;
                S1(sub2ind([tot, ttot], idx, idx - (2 * (boundary_UD < 2*M-1) - 1))) = S1(sub2ind([tot, ttot], idx, idx - (2 * (boundary_UD < 2*M-1) - 1)))+1;
            end

            S1(sub2ind([tot, ttot], idx, idx)) = S1(sub2ind([tot, ttot], idx, idx))-1;
            S1(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1))) = S1(sub2ind([tot, ttot], idx, idx + (2 * (boundary_UD < 2*M-1) - 1)))+1;

        end


        inner_LR = {
            2*M+1:M:tot-3*M+1,                      % edge 1 (e.g., top inner indices)
            3*M:M:tot-2*M            % edge 2 (e.g., bottom inner indices)
            };

        for edge = 1:2
            idx = inner_LR{edge};

            S1(sub2ind([tot, ttot], idx, idx)) = S1(sub2ind([tot, ttot], idx, idx))-2;
            S1(sub2ind([tot, ttot], idx, idx + M)) = S1(sub2ind([tot, ttot], idx, idx + M))+1;
            S1(sub2ind([tot, ttot], idx, idx - M)) = S1(sub2ind([tot, ttot], idx, idx - M))+1;
        end


        inner_UD = {
            3:M-2,                      % edge 1 (e.g., top inner indices)
            tot-N+3 : tot-2            % edge 2 (e.g., bottom inner indices)
            };

        for edge = 1:2
            idx = inner_UD{edge};

            S1(sub2ind([tot, ttot], idx, idx)) = S1(sub2ind([tot, ttot], idx, idx))-2;
            S1(sub2ind([tot, ttot], idx, idx + 1)) = S1(sub2ind([tot, ttot], idx, idx + 1))+1;
            S1(sub2ind([tot, ttot], idx, idx - 1)) = S1(sub2ind([tot, ttot], idx, idx - 1))+1;
        end



        % Collect all boundary rows (top, bottom, left, right)
        top_edge    = 1:M;
        bottom_edge = M*(N-1)+1 : M*N;
        left_edge   = 1:M:N*M;
        right_edge  = M:M:N*M;

        % Combine unique indices
        boundary_rows = unique([top_edge, bottom_edge, left_edge, right_edge]);
        S1(boundary_rows,:) = S1(boundary_rows,:);
        S2 = 0;

    end
  


end

