function [f, g, fl, fr, gl, gr] = ScalingFunction(pm, cont, x0, alfa, sfunparam,N, fl, fr, gl, gr)
    % pm sets the sign of the power of the scaling function
    % cont (0 or 1) determines if continuity is enforced
    % x0 is the location of the discontinuity
    % alfa is the power of the scaling function far from the origin
    % Al and Ar are prefactors
    % s is the difference in shift away from the singularity between left and right
    % N is the number of terms in the exponent
    % v contains the parameters of the functions in the exponent

    % [f,g]=ScalingFunction(-1,1,0.5,1,1,1,1,1,1,[0.5 1 2 1 0.5 1 2 1]);
    % x=linspace(-2,2,100);
    % plot(x,f(x));

    Al=sfunparam(1);
    Ar=sfunparam(2);
    sl=sfunparam(3);
    sr=sfunparam(4);
    v=sfunparam(5:(4 + 8 * N));

    
    if isempty(fl)
        syms x y z t s
        k1 = sym('a', [N, 1]);
        k2 = sym('b', [N, 1]);
        k3 = sym('c', [N, 1]);
        k4 = sym('d', [N, 1]);

        str1 = '(x,y,z,t';
        strl2 = ',s)=abs(y)*abs(x-(z+abs(s)*(pm==-1)))^(pm*abs(t) ';
        strr2 = ',s)=abs(y)*abs(x-(z-abs(s)*(pm==-1)))^(pm*abs(t) ';
        for i = 1:N
            ii = num2str(i);
            str1 = [str1, ',k1(', ii, '),k2(', ii, '),k3(', ii, '),k4(', ii, ')'];
            strl2 = [strl2, '+k1(', ii, ')/N/(abs(x-(z+abs(1+k2(', ii, ')))).^k3(', ii, ') + k4(', ii, '))'];
            strr2 = [strr2, '+k1(', ii, ')/N/(abs(x-(z-abs(1+k2(', ii, ')))).^k3(', ii, ') + k4(', ii, '))'];
        end

        eval(['Fl', str1, strl2, ');']);
        %Fl(x,y,z,t,a1,a2,a3,a4,s)=abs(y)*abs(x-(z+abs(s)*(pm==-1)))^(pm*abs(t) + a1/N/(abs(x-(z+abs(a2))).^a3 + a4) );

        if pm == -1
            eval(['Fr', str1, strr2, ');']);
            %Fr(x,y,z,t,a1,a2,a3,a4,s)=abs(y)*abs(x-(z-abs(s)*(pm==-1)))^(pm*(abs(t) + a1/N/(abs(x-(z-abs(a2))).^a3 + a4) ));
        end

        Gl = diff(Fl, x);

        ffl = matlabFunction(Fl); eval(['fl=@(x,y,z,t,k1,k2,k3,k4,s)ffl', str1, ',s);']);
        ggl = matlabFunction(Gl); eval(['gl=@(x,y,z,t,k1,k2,k3,k4,s)ggl', str1, ',s);']);
        if pm == 1
            fr = @(x, y, z, t, k1, k2, k3, k4, s)0; gr = @(x, y, z, t, k1, k2, k3, k4, s)0;
        else
            Gr = diff(Fr, x);
            ffr = matlabFunction(Fr); eval(['fr=@(x,y,z,t,k1,k2,k3,k4,s)ffr', str1, ',s);']);
            ggr = matlabFunction(Gr); eval(['gr=@(x,y,z,t,k1,k2,k3,k4,s)ggr', str1, ',s);']);
        end
        x = []; y = []; z = []; t = []; s = []; k1 = []; k2 = []; k3 = []; k4 = [];
    end

    alfal = []; betal = []; etal = []; deltal = []; alfar = []; betar = []; etar = []; deltar = [];
    for i = 1:N
        alfal(i) = v(0 + i);
        betal(i) = v(1 + i);
        etal(i) = v(2 + i);
        deltal(i) = v(3 + i);
        alfar(i) = v(4 + i);
        betar(i) = v(5 + i);
        etar(i) = v(6 + i);
        deltar(i) = v(7 + i);
    end

    if cont && pm == -1
        f = @(y)fun(y, x0, @(x)fl(x, Al, x0, alfa, alfal, betal, etal, deltal, sl), @(x)fr(x, Ar, x0, alfa, alfar, betar, etar, deltar, sr) ...
            * fl(x0, Al, x0, alfa, alfal, betal, etal, deltal, sl) / fr(x0, Ar, x0, alfa, alfar, betar, etar, deltar, sr));
        g = @(y)fun(y, x0, @(x)gl(x, Al, x0, alfa, alfal, betal, etal, deltal, sl), @(x)gr(x, Ar, x0, alfa, alfar, betar, etar, deltar, sr) ...
            * fl(x0, Al, x0, alfa, alfal, betal, etal, deltal, sl) / fr(x0, Ar, x0, alfa, alfar, betar, etar, deltar, sr));
    else
        f = @(y)fun(y, x0, @(x)fl(x, Al, x0, alfa, alfal, betal, etal, deltal, sl), @(x)fr(x, Ar, x0, alfa, alfar, betar, etar, deltar, sr));
        g = @(y)fun(y, x0, @(x)gl(x, Al, x0, alfa, alfal, betal, etal, deltal, sl), @(x)gr(x, Ar, x0, alfa, alfar, betar, etar, deltar, sr));
    end

end

function x = fun(y, x0, fl, fr)
    x = y;
    x(y < x0) = fl(y(y < x0));
    x(y >= x0) = fr(y(y >= x0));
end
