function l_cylinder = lenghtTroughCylinder(Q, C_pixel, C_cylinder, R, H)
    N = size(C_pixel,1);
    QCylinder = Q - C_cylinder;
    pixelQ = C_pixel - Q;

    A = sum(pixelQ(:, [1, 3]).^2, 2);
    B = 2*sum((pixelQ(:, [1, 3]) .* QCylinder([1, 3])),2);
    C = sum(QCylinder([1, 3]).^2, 2) - R^2;

    D = B.^2 - 4.*A.*C;

    t1 = nan(N,1);
    t2 = nan(N,1);
    T1 = nan(N,3);
    T2 = nan(N,3);
    l_cylinder = zeros(N,1);
    
    valid = D >= 0 & A > 0;
    sqrtD = sqrt(D(valid));

    t_minus = (-B(valid) - sqrtD) ./ (2 .* A(valid));
    t_plus = (-B(valid) + sqrtD) ./ (2 .* A(valid));

    t1(valid) = min(t_minus, t_plus);
    t2(valid) = max(t_minus, t_plus);

    if size(Q,1)>1
        Q = Q(valid);
    end

    T1(valid, :) = (t1(valid) .* pixelQ(valid, :)) + Q;
    T2(valid, :) = (t2(valid) .* pixelQ(valid, :)) + Q;

    for i = 1:N
        if T1(i, 2) > C_cylinder(2) + H/2
            t1_corr = (H/2 - QCylinder(2))/pixelQ(i, 2);
            T1(i, :) = (t1_corr .* pixelQ(i, :)) + Q;
        elseif T1(i, 2) < C_cylinder(2) - H/2
            t1_corr = (- H/2 - QCylinder(2))/pixelQ(i, 2);
            T1(i, :) = (t1_corr .* pixelQ(i, :)) + Q;
        end
    
        if T2(i, 2) > C_cylinder(2) + H/2
            t2_corr = (H/2 - QCylinder(2))/pixelQ(i, 2);
            T2(i, :) = (t2_corr .* pixelQ(i, :)) + Q;
        elseif T2(i, 2) < C_cylinder(2) - H/2
            t2_corr = (- H/2 - QCylinder(2))/pixelQ(i, 2);
            T2(i, :) = (t2_corr .* pixelQ(i, :)) + Q;
        end
    end

    l_cylinder(valid) = sqrt(sum((T2(valid, :) - T1(valid, :)).^2, 2));

end