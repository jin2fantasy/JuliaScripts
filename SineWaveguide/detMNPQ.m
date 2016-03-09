function result = detMNPQ(f, beta0)
    % constants
    c = 299792458;
    mu0 = 4*pi*10e-7;
    epsilon0 = 1/(mu0*c^2);
    %%%

    % geometry parameters
    a = 0.77e-3;
    p = 0.46e-3;
    h = 0.43e-3;
    g = (0.15e-3)/2;
    s = p/2;
    %%%

    % space harmonics number in region I
    n = -3:3;
    % standing wave number in region II and III
    m = 0:4;
    % number of steps in region II and III
    NN = 1000;
    %%%

    %%%
    omega = 2*pi*f;
    k = omega*sqrt(mu0*epsilon0);
    kx = pi/a; % propagation constant in x direction
    % propagation constant in z direction (region I)
    betan = beta0 + 2*n*pi/p;
    % propagation constant for standing wave in z direction (region II, III)
    % kzmn(m, N)
    kzmn = repmat(m', 1, NN)*pi./repmat((p-p/pi*acos((NN-2*(1:NN)+2)/NN)), numel(m), 1);

    % propagation constant along y direction in region I
    nu = @(f) sqrt(abs((k^2 - betan.^2 - kx^2)));
    nu_index = @(f) ((k^2 - betan.^2 - kx^2)) < 0;
    % propagation constant for standing wave along y direction (region II, III)
    lmn = @(f,nn) sqrt(abs((k^2  - kzmn(:,nn).^2 - kx^2)));
    lmn_index = @(f,nn) ((k^2  - kzmn(:,nn).^2 - kx^2)) < 0;
    % some functions
    Fl = @(f,nn,x) lmn_index(f,nn).*sinh(lmn(f,nn)*x) + ~lmn_index(f,nn).*sind(lmn(f,nn)*x*180/pi);
    Fpl = @(f,nn,x) lmn_index(f,nn).*cosh(lmn(f,nn)*x) + ~lmn_index(f,nn).*cosd(lmn(f,nn)*x*180/pi);
    Gl = @(f,nn,x) lmn_index(f,nn).*cosh(lmn(f,nn)*x) + ~lmn_index(f,nn).*cosd(lmn(f,nn)*x*180/pi);
    Gpl = @(f,nn,x) lmn_index(f,nn).*sinh(lmn(f,nn)*x) + ~lmn_index(f,nn).*(-sind(lmn(f,nn)*x*180/pi));
    Fnu = @(f) nu_index(f).*sinh(nu(f)*g) + ~nu_index(f).*sind(nu(f)*g*180*pi);
    Fpnu = @(f) nu_index(f).*cosh(nu(f)*g) + ~nu_index(f).*cosd(nu(f)*g*180*pi);
    Gnu = @(f) nu_index(f).*cosh(nu(f)*g) + ~nu_index(f).*cosd(nu(f)*g*180*pi);
    Gpnu = @(f) nu_index(f).*sinh(nu(f)*g) + ~nu_index(f).*(-sind(nu(f)*g*180*pi));
    Fl1 =  Fl(f,1,g);
    Fpl1 =  Fpl(f,1,g);
    Gl1 =  Gl(f,1,g);
    Gpl1 =  Gpl(f,1,g);
    lm1 = lmn(f, 1);
    % Rp_index(n, m)
    R_index = abs(repmat(betan', 1, numel(m))) ~= abs(repmat(m*pi/p, numel(n), 1));
    Rp = R_index.*1i.*repmat(betan', 1, numel(m)).*(1-(cosd(repmat(betan', 1, numel(m))*p*180/pi) + 1i*sind(repmat(betan', 1, numel(m))*p*180/pi)).*repmat((-1).^m, numel(n), 1))./((repmat((betan.^2)', 1, numel(m)) - repmat((m*pi/p).^2, numel(n), 1)));
    Rp(~R_index) = p/2;
    Rp(:,1) = Rp(:,1) + ~R_index(:,1)*p/2;
    Rm = R_index.*(-1i).*repmat(betan', 1, numel(m)).*(1-(cosd(repmat(betan', 1, numel(m))*p*180/pi) - 1i*sind(repmat(betan', 1, numel(m))*p*180/pi)).*repmat((-1).^m, numel(n), 1))./((repmat((betan.^2)', 1, numel(m)) - repmat((m*pi/p).^2, numel(n), 1)));
    Rm(~R_index) = p/2;
    Rm(:,1) = Rm(:,1) + ~R_index(:,1)*p/2;

    % calculate bm1/am1 and dm1/cm1
    bmpam = -Fpl(f,NN,g+h)./Gpl(f,NN,g+h);
    dmpcm = Fpl(f,NN,g+h)./Gpl(f,NN,g+h);
    bmpam(isnan(bmpam)) = -1;
    dmpcm(isnan(dmpcm)) = 1;
    for kk = NN:-1:2
        yn1 = g + h*(kk-1)/NN;
        dn = p - p/pi*acos((NN-2*kk+2)/NN);
        dn1 = p - p/pi*acos((NN-2*kk+4)/NN);
        bmpam_index = (bmpam == -1);
        Ymn1n1 = ~bmpam_index.*(Fl(f,kk,yn1) + bmpam.*Gl(f,kk,yn1))./(lmn(f,kk).*Fpl(f,kk,yn1) + bmpam.*Gpl(f,kk,yn1));
        Ymn1n1(isnan(Ymn1n1)) = 0;
        Ymn1n1 = Ymn1n1 - bmpam_index./lmn(f,kk)*dn1/dn;
        bmpam = (Ymn1n1.*lmn(f,kk-1).*Fpl(f,kk-1,yn1) - Fl(f,kk-1,yn1))./(Gl(f,kk-1,yn1) - Ymn1n1.*lmn(f,kk-1).*Gpl(f,kk-1,yn1));
        Ymn1n1 = (Fl(f,kk,yn1) + dmpcm.*Gl(f,kk,yn1))./(lmn(f,kk).*Fpl(f,kk,yn1) + dmpcm.*Gpl(f,kk,yn1))*dn1/dn;
        dmpcm = (Ymn1n1.*lmn(f,kk-1).*Fpl(f,kk-1,yn1) - Fl(f,kk-1,yn1))./(Gl(f,kk-1,yn1) - Ymn1n1.*lmn(f,kk-1).*Gpl(f,kk-1,yn1));
    end
    M = eye(numel(n))*p^2.*repmat(nu(f)', 1, numel(n)).*repmat(Fpnu(f)', 1, numel(n)) - repmat(Fnu(f), numel(n), 1)*lm1(1)*(Fpl1(1) + bmpam(1)*Gpl1(1))/(Fl1(1) + bmpam(1)*Gl1(1)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,1)', numel(n), 1);
    N = eye(numel(n))*p^2.*repmat(nu(f)', 1, numel(n)).*repmat(Gpnu(f)', 1, numel(n)) - repmat(Gnu(f), numel(n), 1)*lm1(1)*(Fpl1(1) + bmpam(1)*Gpl1(1))/(Fl1(1) + bmpam(1)*Gl1(1)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,1)', numel(n), 1);
    P = eye(numel(n))*p^2.*repmat(nu(f)', 1, numel(n)).*repmat(Fpnu(f)', 1, numel(n)) + repmat(Fnu(f), numel(n), 1)*lm1(1)*(Fpl1(1) - dmpcm(1)*Gpl1(1))/(-Fl1(1) + dmpcm(1)*Gl1(1)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,1)', numel(n), 1).*exp(1i*(repmat(betan', 1, numel(n)) - repmat(betan, numel(n), 1))*s);
    Q = -eye(numel(n))*p^2.*repmat(nu(f)', 1, numel(n)).*repmat(Gpnu(f)', 1, numel(n)) - repmat(Gnu(f), numel(n), 1)*lm1(1)*(Fpl1(1) - dmpcm(1)*Gpl1(1))/(-Fl1(1) + dmpcm(1)*Gl1(1)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,1)', numel(n), 1).*exp(1i*(repmat(betan', 1, numel(n)) - repmat(betan, numel(n), 1))*s);
    for kk = 2:numel(m)
        M = M - repmat(Fnu(f), numel(n), 1)*2*lm1(kk)*(Fpl1(kk) + bmpam(kk)*Gpl1(kk))/(Fl1(kk) + bmpam(kk).*Gl1(kk)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,kk)', numel(n), 1);
        N = N - repmat(Gnu(f), numel(n), 1)*2*lm1(kk)*(Fpl1(kk) + bmpam(kk)*Gpl1(kk))/(Fl1(kk) + bmpam(kk).*Gl1(kk)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,kk)', numel(n), 1);
        P = P + repmat(Fnu(f), numel(n), 1)*2*lm1(kk)*(Fpl1(kk) - dmpcm(kk)*Gpl1(kk))/(-Fl1(kk) + dmpcm(kk).*Gl1(kk)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,kk)', numel(n), 1).*exp(1i*(repmat(betan', 1, numel(n)) - repmat(betan, numel(n), 1))*s);
        Q = Q - repmat(Gnu(f), numel(n), 1)*2*lm1(kk)*(Fpl1(kk) - dmpcm(kk)*Gpl1(kk))/(-Fl1(kk) + dmpcm(kk).*Gl1(kk)).*repmat(Rp(:,1), 1, numel(n)).*repmat(Rm(:,kk)', numel(n), 1).*exp(1i*(repmat(betan', 1, numel(n)) - repmat(betan, numel(n), 1))*s);
    end
    result = det([M N; P Q]);
    result = imag(result);
end
