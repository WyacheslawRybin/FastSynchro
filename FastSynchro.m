clear all; close all;

TT = 100;               % Transient time
CT = 20;                % Computation time
WT = 0.15;              % Window time
h = 0.0025;              % Integration step time
a = [0.5 10 28 8/3]';   % parameters, a(1) - symmetry coefficient
Kforward = [0 40 0]';   % Synchronization coefficients for forward synchronization
Kbackward = [0 40 0]';  % Synchronization coefficients for backward synchronization
X = [3 -3 0]';          % Initial conditions for master system
X1 = [1 -3 0]';         % Initial conditions for slave system
itrs = 100;             % Amount of synchronization iterations
y = 20;                 % Final array decimation coefficient

% Transient time calculation
for i = 1:ceil(TT/h)
    X = Lorenz(X,a,h,[0 0 0],[0 0 0]);
    X1 = Lorenz(X1,a,h,[0 0 0],[0 0 0]);
end

X1_start = X1;
Xwrite = zeros(3, ceil(CT/h/y));

% Time domain calculation
m = 0;
for i = 1:ceil(CT/h)
    X = Lorenz(X,a,h,[0 0 0],[0 0 0]);
    if mod(i,y) == 0
        m= m +1;
        Xwrite(:,m)= X;
    end
end

% Initialization of helpfull arrays for calculation
WT_forward = zeros(3, ceil(WT/h));
buffer_norm = zeros(1, ceil(WT/h)-1);
buffer_rms = zeros(1, itrs);
buffer_last_rms = zeros(1, ceil(CT/h/y));
WT_iter = ceil(WT/h);


hw = waitbar(0,'Please wait...');

% Calculation of Forward-Backward synchronization for every y point in time domain
for k = 1:m
    waitbar(k/m,hw,'Processing...');
    X = Xwrite(:,k);
    %X1 = X1_start;
    X1 = X+5 ;
    % Window array calculation
    for i = 1:WT_iter
        WT_forward(:,i) = X;
        X = Lorenz(X,a,h,[0 0 0],[0 0 0]);
    end

    % Formatting window array for backward synchronization
    WT_backward = flip(WT_forward');
    WT_backward = WT_backward';

    %rms_error = log10(err1) - log10(err0);

    for i = 1:itrs
        %Forward synch
        for j = 1:(ceil(WT/h)-1)
            buffer_norm(j) = norm(abs(X1-WT_forward(:,j)'));
            X1 = Lorenz(X1,a,h,WT_forward(:,j),Kforward);
        end
        %Backward synch
        for j = 1:(ceil(WT/h)-1)
            X1 = Lorenz(X1,a,-h,WT_backward(:,j),-Kbackward);
        end

        buffer_rms(i) = rms(buffer_norm);

    end
    buffer_last_rms(k) = log10(buffer_rms(end)) - log10(buffer_rms(1));

end
close(hw);

%isnan checking
buffer_last_rms(isnan(buffer_last_rms)) = 1000000;

figure
surf([Xwrite(1,:)', Xwrite(1,:)'], [Xwrite(2,:)', Xwrite(2,:)'], [Xwrite(3,:)', Xwrite(3,:)'],...
    [buffer_last_rms(1,:)',buffer_last_rms(1,:)'],'EdgeColor','flat', 'FaceColor','none',LineWidth=1.5);
zlabel('$z$','interpreter','latex','FontSize',15);
ylabel('$y$','interpreter','latex','FontSize',15);
xlabel('$x$','interpreter','latex','FontSize',15);
colorbar;
colormap([turbo(1000); 1-flip(copper(144));])
caxis([-14 2]);
view(0,0);

% function Lorenz:
% X - state variables
% a - parameters, a(1) - symmetry coefficient
% h - integration step
% S - master system state variables (for synchronization)
% K - array of synchronization coefficients

function Y = Lorenz(X,a,h,S,K)

%Pecora-Carroll synchronization
N(1) = K(1) * (S(1) - X(1));
N(2) = K(2) * (S(2) - X(2));
N(3) = K(3) * (S(3) - X(3));

h1 = h * a(1); %a(1) - symmetry coefficient
h2 = h * (1 - a(1));

X(1) = X(1) + h1 * ( a(2) * ( X(2) - X(1)) + N(1)) ;
X(2) = X(2) + h1 * ( X(1) * ( a(3) - X(3)) - X(2) + N(2)) ;
X(3) = X(3) + h1 * ( X(1) * X(2) - a(4) * X(3) +N(3)) ;

X(3) = (X(3) + h2 * ( X(1) * X(2) + N(3)))/(1 + h2*a(4)) ;
X(2) = (X(2) + h2 * ( X(1) * ( a(3) - X(3)) + N(2)))/(1 + h2) ;
X(1) = (X(1) + h2 * ( a(2) * ( X(2)) + N(1)))/(1 + a(2)*h2) ;

Y(1) = X(1);
Y(2) = X(2);
Y(3) = X(3);
end