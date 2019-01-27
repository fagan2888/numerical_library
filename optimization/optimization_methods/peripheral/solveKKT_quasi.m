% solve KKT system efficiently
function [pk,y,pn] = solveKKT_quasi(g,B,c,J,opt)

    if strcmp(opt,'tq') == 1
        
        % get system size
        n = numel(g);
        m = numel(c);

        % construct change of variables (JQ = [0 T])
        [Qt,Tt] = qr(J'); % get transpose variables
        P = fliplr(eye(n)); % create permutation matrix
        T = (P*Tt)';T = T(:,n-m+1:end); % get T factor
        Q = Qt*P; % get Q factor
        Z = Q(:,1:n-m); % nullspace of J
        Y = Q(:,n-m+1:end); % rangespace of J'

        % construct QP subproblem solution
        pn = Y*linsolve(T,-c); % feasibility step
        pt = Z*linsolve(Z'*B*Z,-Z'*(g+B*pn)); % minimization step
        pk = pt + pn; % overall solution step

        % construct optimal QP subproblem multiplier
        y = linsolve(T',Y'*(g+B*pk));
        
    elseif strcmp(opt,'pbs') == 1
        
        % get system size
        n = numel(g);
        m = numel(c);
        
        % construct change of variables (JP = [Bj S])
        [~,~,Pt] = lu(J');P = Pt';
        BjS = J*P;
        Bj = BjS(:,1:m);
        S = BjS(:,m+1:end);
        Z = P*[-linsolve(Bj,S);eye(n-m)];
        Y = P*[eye(m);zeros(n-m,m)];
        
        % construct QP subproblem solution
        pn = Y*linsolve(Bj,-c); % feasibility step
        pt = Z*linsolve(Z'*B*Z,-Z'*(g+B*pn)); % minimization step
        pk = pt + pn; % overall solution step
        
        % construct optimal QP subproblem multiplier
        y = linsolve(Bj',Y'*(g+B*pk));
        
    elseif strcmp(opt,'w') == 1
        
        % construct change of variables (JQ = [0 W])
        Z = null(J);
        Y = orth(J');
        W = J*Y;
        
        % construct QP subproblem solution
        pn = Y*linsolve(W,-c); % feasibility step
        pt = Z*linsolve(Z'*B*Z,-Z'*(g+B*pn)); % minimization step
        pk = pt + pn; % overall solution step
        
        % construct optimal QP subproblem multiplier
        y = linsolve(W',Y'*(g+B*pk));
        
    end


end