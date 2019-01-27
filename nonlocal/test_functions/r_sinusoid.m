function f = r_sinusoid(t,p,A,w,operation,type)

    % return appropriate operation
    switch operation
        
        % differentiation
        case 'diff'

            % return appropriate function
            switch type
                case 'sin'
                    scale = A*(t^(1-p))*w;
                    term1 = hypergeom(1,[3/2-p/2,2-p/2],-(1/4)*(t^2)*(w^2))/gamma(2-p);
                    term2 = 2*(t*w)^2*hypergeom(2,[5/2-p/2,3-p/2],-(1/4)*(t^2)*(w^2))/gamma(5-p);
                    f = scale*(term1-term2);
                case 'cos'
                    scale = A*t^(-p);
                    term1 = hypergeom(1,[1-p/2,3/2-p/2],-(1/4)*(t^2)*(w^2))/gamma(1-p);
                    term2 = 2*(t*w)^2*hypergeom(2,[2-p/2,5/2-p/2],-(1/4)*(t^2)*(w^2))/gamma(4-p);
                    f = scale*(term1-term2);
            end
            
        % integration
        case 'int'

            % return appropriate function
            switch type
                case 'sin'
                    f = ((A*(t^(p+1))*w)/gamma(p+2))*hypergeom(1,[1+p/2,3/2+p/2],-(1/4)*(t^2)*(w^2));
                case 'cos'
                    f = ((A*(t^p))/gamma(p+1))*hypergeom(1,[1/2+p/2,1+p/2],-(1/4)*(t^2)*(w^2));
            end
    end

end