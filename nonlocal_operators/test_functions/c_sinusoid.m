function f = c_sinusoid(t,p,A,w,type)

    % return appropriate function
    switch type
        case 'sin'
            f = ((A*w*t^(1-p))/gamma(2-p))*hypergeom(1,[1-p/2,3/2-p/2],-(1/4)*(t^2)*(w^2));
        case 'cos'
            f = -((A*(w^2)*t^(2-p))/gamma(3-p))*hypergeom(1,[3/2-p/2,2-p/2],-(1/4)*(t^2)*(w^2));
        case 'cos2'
            f = -((A*t^(2-p)*w^2)/gamma(3-p))*hypergeom(1,[3/2-p/2,2-p/2],-(1/4)*(t^2)*(w^2));
    end

end