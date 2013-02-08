function inf = playiterateFB(model)
% function inf = playiterateFB(model)

L = model.seqlength;
p = model.bias;
e = model.err;
N = model.nreads;
dna = model.dna;

reads = model.reads;

% initialize S estimate matrix
S = 0.25* ones(L,4);
S(1,:) = [0 0 0 0];
S(2,:) = [0 0 0 0];
S(end,:) = [0 0 0 0];
S(end-1,:) = [0 0 0 0];
S(1,dna(1)) = 1;
S(2,dna(2)) = 1;
S(end,dna(end)) = 1;
S(end-1,dna(end-1)) = 1;

SHold = S;
% build transition matrix
A = gen_transmatrix(L,p);

for i=1:N
    T(i) = length(reads(i).z);
    alpha{i} = zeros(T(i),L);
    alpha{i}(1,1) = 1;
    alpha{i}(2,2) = 1;
    beta{i} = zeros(T(i),L);
end

qcnt = 1;
maxiter = 10;
epsilon = 0.001;

count = 1;

q = 2; % EM iteration

for iter = 2:max(T)
    clear sllhd;
    sllhd(1) = -Inf;
    % do EM for each read up to time _iter_
    
    while q < maxiter
        fprintf('doing em iteration %d at t=%d\n',q-1,iter);
        for i=1:N
            if iter<=T(i)
                % E-step for read i
                
                %alpha
                for t = 2:iter
                    alpha{i}(t,:) = alpha{i}(t-1,:) * A * diag(S(:,reads(i).x(t)));
                    c{i}(t) = sum(alpha{i}(t,:));
                    alpha{i}(t,:) = alpha{i}(t,:) / c{i}(t);
                end

                %beta
                beta{i}(iter,:) = ones(1,L);
                for t=iter-1:-1:1
                    beta{i}(t,:) = A * diag(S(:,reads(i).x(t+1))) * beta{i}(t+1,:)';
                    beta{i}(t,:) = beta{i}(t,:) / c{i}(t+1);
                end

                % gamma
                gamma{i} = alpha{i} .* beta{i};
                gamma{i} = gamma{i} / sum(gamma{i}(1,:));
                gamma{i}(isnan(gamma{i})==1) = 0;

                lc = log(c{i});
                lc(lc==-Inf) = 0;
                lc(isnan(lc)==1) = 0;
                llhd(i) = sum(lc);

                % M-step

                tmp = sum(gamma{i},1);
                SS{i} = S; % just controls for where you haven't visited yet
                for d = 1:4
                    SS{i}(:,d) = (sum(gamma{i}(model.reads(i).x==d,:),1) ./ tmp)';
                end
                SS{i}(isnan(SS{i}==1)) = 0.25;
                SS{i} = SS{i} ./ repmat(sum(SS{i},2),1,4); % normalize
            end
        end

        % average S
        S = mean(cat(3, SS{:}), 3);
        tmp = isnan(S);
        S(isnan(S)) = 0;
        S = S + tmp.*SHold;
        inf.h(count).S = S;
        inf.h(count).inf_ent = calc_entropy(S);
        count = count+1;
        %pause

        % convergence check
        
        sllhd(q) = sum(llhd);

        fprintf('t =%d, em-iter = %d, sllhd = %e\n',iter,q,sllhd(q));
        fprintf('sllhd - sllhdprev = %e\n',abs(sllhd(q) - sllhd(q-1))); 
        if abs(sllhd(q) - sllhd(q-1)) < epsilon*abs(sllhd(q)) | sllhd(q)==0
            break
        end

        q = q + 1;
        %pause 
    end

    q = 2; % reset q
    clear llhd
end

inf.sllhd = sllhd;
inf.alpha = alpha;
inf.beta = beta;
inf.c = c;
inf.gamma = gamma;
inf.SS = SS;
inf.S = S;
inf.A = A;
