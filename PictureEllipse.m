% Parameters
N     = 10000;       % total length
p     = 0.4;         % proportion of 1's
imax  = 50;          % number of simulations to plot
h     = 3.841;       % chi-square(1) threshold

nOnes  = round(p*N);
nZeros = N - nOnes;

figure('Color','w'); hold on;

% Plot multiple random walks (semi-transparent)
alphaVal = 0.15;
for sim = 1:imax
    x = zeros(1,N+1);
    y = zeros(1,N+1);
    remainingZeros = nZeros;
    remainingOnes  = nOnes;
    for k = 1:N
        px = remainingZeros / (remainingZeros + remainingOnes);
        if rand < px
            x(k+1) = x(k) + 1;
            y(k+1) = y(k);
            remainingZeros = remainingZeros - 1;
        else
            x(k+1) = x(k);
            y(k+1) = y(k) + 1;
            remainingOnes = remainingOnes - 1;
        end
    end
    % Plot with low alpha by plotting a faded RGB color
    plot(x, y, 'Color',[0 0 1 alphaVal], 'LineWidth', 1);
end

% Analytic G=0 curve via solving quadratic in s for each x
A = nZeros;
B = nOnes;

% Coefficients that do not depend on x:
alpha = A^2 + (h*A*B)/N;
const_hAB = h*A*B;  % reuse

% Choose resolution for x sweep (fine enough to look smooth)
xvec = linspace(0, nZeros, 2000);  

s1 = nan(size(xvec));
s2 = nan(size(xvec));

for ii = 1:numel(xvec)
    xval = xvec(ii);
    beta = - (2*A*N*xval + const_hAB);
    gamma = (N^2) * (xval^2);
    disc = beta^2 - 4*alpha*gamma;
    if disc < 0
        continue;  % no real solution at this x
    end
    sqrtD = sqrt(disc);
    sSol1 = (-beta + sqrtD) / (2*alpha);
    sSol2 = (-beta - sqrtD) / (2*alpha);
    % Keep only solutions in [0, N]
    if sSol1 >= 0 && sSol1 <= N
        s1(ii) = sSol1;
    end
    if sSol2 >= 0 && sSol2 <= N
        s2(ii) = sSol2;
    end
end

% Compute y from s - x and keep only physically valid points:
y1 = s1 - xvec;
y2 = s2 - xvec;

valid1 = ~isnan(y1) & (y1 >= 0) & (y1 <= nOnes);
valid2 = ~isnan(y2) & (y2 >= 0) & (y2 <= nOnes);

% Plot branches
plot(xvec(valid1), y1(valid1), 'r-', 'LineWidth', 2);
plot(xvec(valid2), y2(valid2), 'r-', 'LineWidth', 2);

xlabel('Cumulative zeros (x)');
ylabel('Cumulative ones (y)');
title(sprintf('Random walks (%d runs) with analytic G=0 curve', imax));
axis equal;
xlim([0 nZeros]);
ylim([0 nOnes]);
grid on;
box on;
legend({'random walks','G=0 curve'}, 'Location','best');
