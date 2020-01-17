close all
clear all
clc

%% input from clustering datasets http://cs.joensuu.fi/sipu/datasets/

% uploads complex vectors
v = importdata('jain.txt')*[1;1i;0]; 

% show them
figure
plot(v,'x')
grid
title('A.K. Jain''s toy problem')

% extract adjacency matrix (distances)
A = 1*(abs(v-v.')<=4); % activate only if distance <=4 !
A = A-diag(diag(A));
N = size(A,1);


%% spectral approach (taken from Lab5, and made compact)

% Laplacians
d = full(sum(A)); % degree vector
D = sum(d); % degrees sum
I = speye(N); % identity matrix
Di = diag(sparse(sqrt(1./d))); % diagonal degrees square-rooted
L = I - Di*A*Di; % normalized Laplacian
M = A*Di*Di; % normalized adjacency matrix
% extract Fiedler`s eigenvectors
[V,DD] = eigs(L,2,'SA');
v1 = Di*V(:,2);
v1 = v1/norm(v1);
% sort
[v1s,pos] = sort(v1,'descend');
% extract the conductance measure by a sweep
conduct = conductance_sweep(A(pos,pos)); 
% identify the minimum -> threshold
[~,mpos] = min(conduct);
threshold = mean(v1s(mpos:mpos+1));
% communities
C1s = sort(pos(1:mpos));
C2s = sort(pos(mpos+1:end));

% show the result
figure
plot(v(C1s),'x','Color',[0 75 135]/255)
hold on
plot(v(C2s),'x','Color',[255,140,0]/255)
hold off
grid
title('Spectral clustering')

disp('spectral approach')
disp(['   Minimum conductance: ' num2str(conduct(mpos))])
disp(['   Cheeger''s lower bound: ' num2str(.5*DD(2,2))])
disp(['   Cheeger''s upper bound: ' num2str(sqrt(2*DD(2,2)))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Community size #1: ' num2str(mpos)])
disp(['   Community size #2: ' num2str(N-mpos)])
disp([' '])


%% PageRank-nibble approach (infinite precision)

i = 1;

% compute the rankink values r1 with teleport set S={i}
% by exploiting Matlab's linear system solution A\b
c = 0.85;
q = zeros(N,1);
q(i) = 1; % teleport vector
r1 = (I-c*M)\((1-c)*q); % ranking vector

% sweep wrt the ordering identified by r1
[u1s,pos2] = sort(r1,'descend');
Nmax = find(u1s>0,1,'last'); % discard nodes with 0 values (never used in Push)
conduct = conductance_sweep(A(pos2,pos2)); 
conduct = conduct(1:Nmax-1); 

% plot the conductance sweep
figure
plot(conduct)
grid
title('PageRank-nibble approach (exact)')

% identify the minimum -> threshold
[~,mpos2] = min(conduct);
threshold2 = mean(u1s(mpos2:mpos2+1));

% identify communities
C1 = sort(pos2(1:mpos2));
C2 = sort(pos2(mpos2+1:end));

% show the result
figure
plot(v(C1),'x','Color',[0 75 135]/255)
hold on
plot(v(C2),'x','Color',[255,140,0]/255)
plot(v(i),'xr')
hold off
grid
title('PageRank-nibble approach (exact)')

disp('PageRank-nibble approach (exact)')
disp(['   Minimum conductance: ' num2str(conduct(mpos2))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Community size #1: ' num2str(mpos2)])
disp(['   Community size #2: ' num2str(N-mpos2)])
disp([' '])


%% PageRank-nibble approach (infinite precision)

i = 1;
epsilon = 1e-3;

% compute the rankink values r2 with teleport set S={i}
% by exploiting the approximate pagerank algorithm with 
% precision epsilon
[r2,complexity] = approx_pagerank(A,q,epsilon,c);

% sweep wrt the ordering identified by r2
[u2s,pos2] = sort(r2,'descend');
Nmax = find(u2s>0,1,'last'); % discard nodes with 0 values (never used in Push)
conduct = conductance_sweep(A(pos2,pos2)); 
conduct = conduct(1:Nmax-1); 

% plot the conductance sweep
figure
plot(conduct)
grid
title('PageRank-nibble approach (approx)')

% identify the minimum -> threshold
[~,mpos2] = min(conduct);
threshold2 = mean(u2s(mpos2:mpos2+1));

% identify communities
C1 = sort(pos2(1:mpos2));
C2 = sort(pos2(mpos2+1:end));


% show the resulting communities
figure
plot(v(C1),'x','Color',[0 75 135]/255)
hold on
plot(v(C2),'x','Color',[255,140,0]/255)
plot(v(i),'xr')
hold off
grid
title('PageRank-nibble approach (approx)')

disp('PageRank-nibble approach (approx)')
disp(['   complexity/D: ' num2str((complexity/D))])
disp(['   epsilon: ' num2str(epsilon)])
disp(['   prec: ' num2str(norm(r1-r2,1))])
disp(['   Minimum conductance: ' num2str(conduct(mpos2))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Community size #1: ' num2str(mpos2)])
disp(['   Community size #2: ' num2str(N-mpos2)])



%% functions

function conduct = conductance_sweep(A)

a = sum(triu(A));
b = sum(tril(A));
D = sum(a+b);
assoc = cumsum(a+b);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:end-1); 

end

function [u, complexity] = approx_pagerank(A,q,ep,c)

N = length(q);
d = sum(A).';
M = A./d'; % <-- corrected
D = sum(d);
u = zeros(N,1); % starting point
v = q; % starting point
th = full(ep*d/D)'; % thresholds
count = 0; % exit counter
complexity = 0; % complexity value (# of operations)
ii = find(q==1); % starting index used for Push operation
while (count<N)
    if v(ii)>th(ii) % push if above threshold
        tmp = v(ii);
        u(ii) = u(ii)+(1-c)*tmp;
        v(ii) = 0;
        v = v + c*M(:,ii)*tmp;    
        complexity = complexity + d(ii); % update complexity
        count = 0; % reset the exit counter
    else % go to next entry if below threshold
        count = count + 1; % increase the exit counter
        ii = mod(ii,N)+1; % update the index used for Push
    end
end

end
