function l_sphere = lenghtTroughSphere(Q, C_pixel, C_sphere, R)

    N = size(C_pixel,1);
    pixelQ = C_pixel - Q;
    QSphere = Q - C_sphere;
    
    A = sum(pixelQ.^2, 2);
    B = 2*sum((pixelQ .* QSphere),2);
    C = sum(QSphere.^2, 2) - R^2;
    
    D = B.^2 - 4.*A.*C;
    
    t1 = nan(N,1);
    t2 = nan(N,1);
    T1 = nan(N,3);
    T2 = nan(N,3);
    l_sphere = zeros(N,1);
    
    valid = D >= 0 & A > 0;
    sqrtD = sqrt(D(valid));
    
    t_minus = (-B(valid) - sqrtD) ./ (2 .* A(valid));
    t_plus = (-B(valid) + sqrtD) ./ (2 .* A(valid));
    
    t1(valid) = min(t_minus, t_plus);
    t2(valid) = max(t_minus, t_plus);

    if size(Q,1)>1
        Q = Q(valid);
    end
    
    T1(valid, :) = t1(valid) .* pixelQ(valid, :) + Q;
    T2(valid, :) = t2(valid) .* pixelQ(valid, :) + Q;
    
    l_sphere(valid) = sqrt(sum((T2(valid, :) - T1(valid, :)).^2, 2));

end