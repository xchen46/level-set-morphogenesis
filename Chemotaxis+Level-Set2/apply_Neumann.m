function [S1, S2] = apply_Neumann(S, M, N, mode)
% APPLY_NEUMANN  Apply Neumann boundary condition patches to sparse matrix
%
% Inputs:
%   S     - Initial sparse operator matrix (M*N x M*N)
%   M, N  - Grid size (rows, columns)
%   mode  - 'symmetric', 'one_sided', or fallback custom mode
%
% Outputs:
%   S1    - Operator patched along left/right boundaries
%   S2    - Operator patched along top/bottom boundaries

    if nargin < 4
        mode = 'symmetric';
    end

    [tot, ttot] = size(S);
    S1 = S;
    S2 = S;

    %% --- Mode 1: One-sided Neumann patch ---
    if strcmp(mode, 'one_sided')
        % Left and Right boundaries
        for boundary_LR = [M+1, tot - 2*M + 1]
            idx = boundary_LR + 1 : boundary_LR + M - 2;

            direction = (2 * (boundary_LR <= M+1) - 1);  % -1 for left, +1 for right
            S1(sub2ind([tot, ttot], idx, idx)) = S1(sub2ind([tot, ttot], idx, idx)) - 1;
            S1(sub2ind([tot, ttot], idx, idx + direction * M)) = S1(sub2ind([tot, ttot], idx, idx + direction * M)) + 1;
        end

        % Top and Bottom boundaries
        for boundary_UD = [2*M+1, 3*M-2]
            idx = boundary_UD + 1 : M : tot - 2*M;
            direction = (2 * (boundary_UD <= 2*M+1) - 1);  % -1 for top, +1 for bottom
            S2(sub2ind([tot, ttot], idx, idx)) = S2(sub2ind([tot, ttot], idx, idx)) - 1;
            S2(sub2ind([tot, ttot], idx, idx + direction)) = S2(sub2ind([tot, ttot], idx, idx + direction)) + 1;
        end
    end

    %% --- Mode 2: Symmetric Neumann patch ---
    if strcmp(mode, 'symmetric')
        % Left and Right boundaries
        for boundary_LR = [M+1, tot - 2*M + 1]
            idx = boundary_LR + 1 : boundary_LR + M - 2;
            direction = (2 * (boundary_LR <= M+1) - 1);  % -1 for left, +1 for right

            k  = sub2ind([tot, ttot], idx, idx);
            kk = sub2ind([tot, ttot], idx, idx + direction * M);
            S1(k)  = S1(k) - 1;
            S1(kk) = S1(kk) + 1;
        end

        % Top and Bottom boundaries
        for boundary_UD = [M+2, 2*M-1]
            idx = boundary_UD : M : tot - M;
            direction = (2 * (boundary_UD < 2*M-1) - 1);  % +1 for top, -1 for bottom

            S2(sub2ind([tot, ttot], idx, idx)) = S2(sub2ind([tot, ttot], idx, idx)) - 1;
            S2(sub2ind([tot, ttot], idx, idx + direction)) = S2(sub2ind([tot, ttot], idx, idx + direction)) + 1;
        end
    end

    %% --- Mode 3: Fallback patch for custom or unknown mode ---
    valid_modes = {'symmetric', 'one_sided'};
    if ~ismember(mode, valid_modes)
        % Left and Right boundaries
        for boundary_LR = [M+1, tot - 2*M + 1]
            idx = boundary_LR : boundary_LR + M - 1;
            direction = (2 * (boundary_LR <= M+1) - 1);  % -1 for left, +1 for right

            k  = sub2ind([tot, ttot], idx, idx);
            kk = sub2ind([tot, ttot], idx, idx + direction * M);

            if boundary_LR == M+1
                S1(k)  = S1(k) + 1;
                S1(kk) = S1(kk) - 1;
            else
                S1(k)  = S1(k) - 1;
                S1(kk) = S1(kk) + 1;
            end
        end

        % Top and Bottom boundaries
        for boundary_UD = [2, M-1]
            idx = boundary_UD : M : tot;

            if boundary_UD == M+2
                S2(sub2ind([tot, ttot], idx, idx))     = S2(sub2ind([tot, ttot], idx, idx)) - 1;
                S2(sub2ind([tot, ttot], idx, idx + 1)) = S2(sub2ind([tot, ttot], idx, idx + 1)) + 1;
            else
                S2(sub2ind([tot, ttot], idx, idx))     = S2(sub2ind([tot, ttot], idx, idx)) + 1;
                S2(sub2ind([tot, ttot], idx, idx - 1)) = S2(sub2ind([tot, ttot], idx, idx - 1)) - 1;
            end
        end
    end

end