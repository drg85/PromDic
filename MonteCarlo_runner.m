% Parameters
N = input('Enter N value: ')
p = input('Enter p value: ')      % proportion of 1's
imax  = 1000;        % number of simulations
h     = 3.841;       % chi-square(1) threshold

nOnes  = round(p*N);
nZeros = N - nOnes;

% Track whether each run has a FP
hasFP    = false(1,imax);   % any G > 0
hasFPEx  = false(1,imax);   % any G > 0 & all counts >= 5

for sim = 1:imax
    %--- Generate sequential random vector (biased walk)
    x = zeros(1,N+1);   % cumulative zeros
    y = zeros(1,N+1);   % cumulative ones
    remainingZeros = nZeros;
    remainingOnes  = nOnes;

    for k = 1:N
        px = remainingZeros / (remainingZeros + remainingOnes);
        if rand < px
            % place a zero
            x(k+1) = x(k) + 1;
            y(k+1) = y(k);
            remainingZeros = remainingZeros - 1;
        else
            % place a one
            x(k+1) = x(k);
            y(k+1) = y(k) + 1;
            remainingOnes = remainingOnes - 1;
        end
    end

    %--- Evaluate G-statistic at all split points
    for k = 1:N-1
        a = y(k+1);        % ones on the left
        b = x(k+1);        % zeros on the left
        c = nOnes - a;     % ones on the right
        d = nZeros - b;    % zeros on the right

        g = (N * (a*d - b*c)^2) / ((a+b)*(c+d)*(a+c)*(b+d)) - h;
        minv = min([a b c d]);

        if g > 0
            hasFP(sim) = true;
            if minv >= 5
                hasFPEx(sim) = true;
            end
        end
    end
end

%--- Results
totalFP   = sum(hasFP);
totalFPEx = sum(hasFPEx);

fprintf('Out of %d runs:\n', imax);
fprintf(' (a) Runs with any G>0: %d (%.2f%%)\n', ...
        totalFP,   100*totalFP/imax);
fprintf(' (b) Runs with any G>0 & all counts>=5: %d (%.2f%%)\n', ...
        totalFPEx, 100*totalFPEx/imax);
