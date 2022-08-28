function [norm_val] = normalize(val, l_min, l_max)
    norm_val = (val - l_min)/((l_max - l_min) + %eps);
    if norm_val > 1 then
        norm_val = 1;
    elseif norm_val < 0 then
        norm_val = 0;
    end
endfunction

function [val]=denormalize(norm_val, l_min, l_max)
    val = norm_val * (l_max - l_min) + l_min;
    if val > l_max then
        val = l_max;
    elseif val < l_min then
        val = l_min;
    end
endfunction

function [dxdt]=dynamicModel(t,c)
    // c(1,:) x position
    //c(2,:) y position
    //c(3,:) movement
    //c(4,:) orientation
    //c(5,:) speed
    //c(6,:) angular speed

    //Parameters
    m = .38;    //mass
    Ip = 0.005; //inertia moment
    d = 0.02;   //distance centroide to axis wheel
    r = 0.03;   //radio wheel
    R = 0.05;   //distance wheel-center

    M = [m 0 ; 0 (Ip+m*d^2)];          //Inertial matrix
    H = [-m*d*c(6)^2 ; m*d*c(5)*c(6)]; //Coriolis and centrifugal forces
    B = [1/r 1/r ; R/r -R/r];          //Conversion  torque-wheel-movil force
    A = [r/2 r/2 ; r/(2*R) -r/(2*R)];  //Speed ratio
    F = [0; 0];                        //Friction
    Ts = [.434 0 ; 0 .434];
    Ks = [2.745 0 ; 0 2.745];
    Kl = [1460.2705 0 ; 0 1460.2705];

    dxdt = [[cos(c(4)) -d*sin(c(4)); sin(c(4)) d*cos(c(4))]*[c(5); c(6)]; c(5); c(6); inv(M + B*inv(Kl)*Ts*inv(A)) * (B*inv(Kl)*Ks*u - (H + F + B*inv(Kl)*inv(A)*[c(5); c(6)]))];
endfunction

function [c,t]= mov(u,ci)
    ti = 0;              //Initial time
    tf = 1;              //Final time
    tspan = [ti:0.1:tf]; //Time vector
    c = ode(ci,ti,tspan,dynamicModel);
    t = tspan;
endfunction

function [ReportPreys, ReportPredators]= ppBehavior(Predators, Preys, dRd, dOd, dAd)
    res = zeros(3,1);                       //Center of mass and dispersion report
    Tv = Predators + Preys;                 //Total of vehicles
    C = zeros(Tv, 6);
    for v=1:Tv
        C(v,4) = denormalize(rand(1,'uniform'),0,2*%pi); ;
    end                                     //Vehicle states
    iPreys = zeros(Preys,3)                 //Influence by prey
    uO=0; uR =0; uA = 0; uI =0;             //Initial speeds (voltages)
    ud = zeros(Tv,2)+uO;                    //Desired speed (voltage)
    lims= 10;                               //Test area limits
    dInf = 4;                               //Maximum distance of prey influence 2
    cp = 0;                                 //Collisions by prey
    m=0;                                    //Time
    collisions = 3;
    dirExp = C(:,4);
    temp = zeros(Tv,1)
    count = zeros(Tv,1)

    for i=1: Preys
        iPreys(i,1) = C(i,1);
        iPreys(i,2) = C(i,2);
        iPreys(i,3) = 0;
    end

    dR = 0.075+dRd; dRp = 0.25; //Respulsion ratio
    dO = 0.075+dOd; dOp = 0.25; //Orientation ratio
    dA = 0.075+dAd; dAp = 0.25; //Atraction ratio

    for v = 1:Tv  //Initialize predators
        if v == 1 then
            C(v,1) = lims*denormalize(rand(1,'uniform'),0.4,0.6); //x position
            C(v,2) = lims*denormalize(rand(1,'uniform'),0.5,0.6); //y position
        else
            while %T
                if v <= Preys then
                    safe = 0;
                    C(v,1) = lims*denormalize(rand(1,'uniform'),0.4,0.6);     //x position
                        C(v,2) = lims*denormalize(rand(1,'uniform'),0.5,0.6); //y position
                    for j = 1:v-1
                        dother = sqrt((C(v,1)-C(j,1))^2 + (C(v,2)-C(j,2))^2);
                        if dother > 0.3 then
                                safe = safe + 1;
                        end
                    end
                    if safe == j then
                        break;
                    end
                else
                    safe = 0;
                    C(v,1) = lims*denormalize(rand(1,'uniform'),0.33,0.66); //x position
                    C(v,2) = lims*denormalize(rand(1,'uniform'),0,0.33);    //y position
                    for j = 1:v-1
                        dother = sqrt((C(v,1)-C(j,1))^2 + (C(v,2)-C(j,2))^2);
                        if dother > 0.3 then
                            safe = safe + 1;
                        end
                    end
                    if safe == j then
                        break;
                    end
                end
            end
        end
        C(v,3) = 0;                                      //Movement made
        C(v,4) = denormalize(rand(1,'uniform'),0,2*%pi); //Orientation
        C(v,5) = 0;                                      //Speed
        C(v,6) = 0;                                      //Angular speed
    end

    scf();
    f = get("current_figure");
    f.BackgroundColor = [1 1 1];
    f.figure_position = [10,10];
    f.figure_size=[700,800];

    while cp < collisions
        m = m + 1;
        for v=1:Tv
            if v <= Preys then
                dR = dRp; dO = dOp; dA = dAp;
            else
                dR = dRd; dO = dOd; dA = dAd;
            end

            //Elements detected
            DetR = 0;      //Repulsion
            DetO = 0;      //Orientation
            DetA = 0;      //Attraction
            DetI = 0;      //Influence (Objects)
            ElemRx = [];   //Elements in repulsion
            ElemOx = [];   //Elements in orientation
            ElemAx = [];   //Elements in attraction
            ElemIx = [];   //Elements in influence
            ElemRy = [];   //Elements in repulsion
            ElemOy = [];   //Elements in orientation
            ElemAy = [];   //Elements in attraction
            ElemIy = [];   //Elements in influence

            //Detection range for repulsion, orientation, attraction, influence, nest and objects
            rangR = 2.0944;
            rangO = 2.0944;
            rangA = 0.5235988;
            rangI = 2.44346;

            //Verify each sensor for repulsion of walls
            for i = 1:5
                if i == 1 then
                    dirObs = C(v,4) -3.8397244;
                else
                    dirObs = dirObs+1.9198622;
                end
                if dirObs < 0 then
                    dirObs = dirObs+(2*%pi);
                elseif dirObs > (2*%pi) then
                    dirObs = dirObs-(2*%pi);
                end
                Dir = [cos(dirObs),sin(dirObs)];
                limitX = C(v,1) + (Dir(1)*dR);
                limitY = C(v,2) + (Dir(2)*dR);

                //Resulting direction due exploration
                if limitX > lims-(dR+0.1) | limitX < (0+dR+0.1) | limitY > lims-(dR+0.1) | limitY < (0+dR+0.1) then
                    explore(v) = 0;
                    dirExp(v) = dirObs+(3*%pi/4)+(rand()*%pi/2);
                    if dirExp(v) > (2*%pi) then
                        dirExp(v) = dirExp(v)-(2*%pi);
                    end
                end
            end

            for z=1:Tv
                if v <> z then //It must not be the same
                    //Angle of the individual with respect to other members of the swarm
                    ang = atan(C(z,2)-C(v,2),C(z,1)-C(v,1))
                    if ang<0 then
                        ang=ang+(2*%pi);
                    elseif ang > (2*%pi) then
                        ang=ang-(2*%pi);
                    end

                    //Calculation of angles of repulsion and attraction with respect to other individuals
                    Beta = ang - C(v,4);
                    if  Beta < 0 then
                        Beta = Beta + (2*%pi);
                    end

                    Gama = C(v,4) - ang;
                    if  Gama < 0 then
                        Gama = Gama + (2*%pi);
                    end

                    if Gama<Beta then
                        Delta = Gama;
                    else
                        Delta = Beta;
                    end


                    //Calculation of the repulsion distance with respect to other individuals
                    if Delta < rangR/2 then
                        dr = sqrt((C(v,1)-C(z,1))^2+(C(v,2)-C(z,2))^2);
                    else
                        dr = %inf;
                    end

                    //Calculation of the orientation distance with respect to objects
                    if Delta < rangO/2 then
                        dor = sqrt((C(v,1)-C(z,1))^2+(C(v,2)-C(z,2))^2);
                    else
                        dor = %inf;
                    end

                     //Calculation of the attraction distance with respect to other individuals
                    if Delta < rangA/2 then
                        da = sqrt((C(v,1)-C(z,1))^2+(C(v,2)-C(z,2))^2);
                    else
                        da = %inf;
                    end

                    //Number of individuals detected in the radius of repulsion, orientation and attraction
                    if  dr <= dR then
                        ElemRx = [ElemRx; cos(ang)];
                        ElemRy = [ElemRy; sin(ang)];
                        DetR = DetR+1; //A vehicle is detected (repuslsion)
                    end

                    if  dor > dR & dor <= dO & DetR == 0 & DetA == 0 then
                        ElemOx = [ElemOx; cos(C(z,4))];
                        ElemOy = [ElemOy; sin(C(z,4))];
                        DetO = DetO+1; //A vehicle is detected (orientation)
                    end

                    if  da > dO & da <= dA & DetR == 0 then
                        ElemAx = [ElemAx; cos(ang)];
                        ElemAy = [ElemAy; sin(ang)];
                        DetA = DetA+1; //A vehicle is detected (attraction)
                    end
                end
            end //of Z

            //Influence prey angle
            if v > Preys then
                for p=1 : Preys
                    if iPreys(p,3) < collisions then
                        dirP(p) = atan(iPreys(p,2)-C(v,2),iPreys(p,1)-C(v,1));
                        if dirP(p)<0 then
                            dirP(p)=dirP(p)+(2*%pi);
                        elseif dirP(p) > (2*%pi) then
                            dirP(p)=dirP(p)-(2*%pi);
                        end

                        //Calculation of influence angle
                        BetaI = dirP(p) - C(v,4);
                            if  BetaI < 0 then
                                BetaI = BetaI + (2*%pi);
                            end

                        GamaI = C(v,4) - dirP(p);
                        if  GamaI < 0 then
                            GamaI = GamaI + (2*%pi);
                        end

                        if GamaI<BetaI then
                            DeltaI = GamaI;
                        else
                            DeltaI = BetaI;
                        end

                        //Calculation of the influence distance
                        if DeltaI<rangI then
                            di(p)=sqrt((iPreys(p,1)-C(v,1))^2+(iPreys(p,2)-C(v,2))^2);
                        else
                            di(p) = %inf;
                        end
                    else
                        di(p) = %inf;
                    end
                end

                minPrey = find(di==min(di))
                if  di(minPrey) <= dInf then
                    ElemIx = [ElemIx; cos(dirP(minPrey))];
                    ElemIy = [ElemIy; sin(dirP(minPrey))];
                    DetI = DetI+1; //A vehicle is detected (repuslsion)
                end

                Ilim = 0.05;
                if (di(minPrey)>=Ilim) then
                    di(minPrey) = di(minPrey);
                else
                    di(minPrey) = Ilim;
                end

                [norm_val]=normalize(di(minPrey),Ilim,dInf)
                [val]=denormalize(norm_val,1.2,2)
                uI = val;
            end

            uO=uI; uR=1; uA=(uI+2)/2;
            uDO=1.2; uDR=1.2; uDA=(uDO+1.2)/2;

            //Average of detected elements
            dirR = atan(-sum(ElemRy),-sum(ElemRx));
            if dirR < 0 then
                dirR = dirR+(2*%pi);
            end

            dirO = atan(sum(ElemOy),sum(ElemOx));
            if dirO < 0 then
                dirO = dirO+(2*%pi);
            end

            dirA = atan(sum(ElemAy),sum(ElemAx));
            if dirA < 0 then
                dirA = dirA+(2*%pi);
            end

            dirI = atan(sum(ElemIy),sum(ElemIx));
            if dirI < 0 then
                dirI = dirI+(2*%pi);
            end

            //Behavior Policies
            //Repulsion rules
            if DetR > 0 then
                if (v > Preys) then
                    ud(v,:) = zeros(1,2)+uR;
                    xT = cos(dirR);
                    yT = sin(dirR);
                    C(v,4) = atan(yT,xT);
                    dirExp(v) = dirR;
                else
                    if iPreys(v,3) < collisions
                        if count(v) < 6 & m > 10 then
                            temp(v) =  temp(v) + 1;
                        else
                            temp(v) = 0;
                            count(v) = 0;
                        end
                        ud(v,:) = zeros(1,2)+uR;
                        xT = cos(dirR);
                        yT = sin(dirR);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirR;
                    end
                end
            end
            //Collisions by prey
            if temp(v) >= 1 then
                count(v) = count(v) + 1;
            end
            if temp(v) == collisions then
                iPreys(v,3) = collisions;
            end
            cp =  mean(iPreys(:,3));

            //Influence rules
            if v > Preys then
                if (DetI>0 & DetR==0) then
                    if DetA == 0 then
                        ud(v,:) = zeros(1,2)+uO;
                        xT = cos(dirI);
                        yT = sin(dirI);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirI;
                    else
                        ud(v,:) = zeros(1,2)+uO;
                        xT = cos(dirI)*0.5 + cos(dirA);
                        yT = sin(dirI)*0.5 + sin(dirA);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirI;
                    end
                end
            end

            //Atraction rules
            if (DetA>0 & DetR==0 & DetO==0 & DetI==0 ) then
                if (v > Preys) then
                    ud(v,:) = zeros(1,2)+uA;
                    xT = cos(dirA);
                    yT = sin(dirA);
                    C(v,4) = atan(yT,xT);
                    dirExp(v) = dirA;
                else
                    if iPreys(v,3) < collisions then
                        ud(v,:) = zeros(1,2)+uDA;
                        xT = cos(dirA);
                        yT = sin(dirA);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirA;
                    end
                end
            end

            //Orientation rules
            if (DetO>0 & DetR==0 & DetA==0 & DetI==0) then
                if (v > Preys) then
                    ud(v,:) = zeros(1,2)+uO;
                    xT = cos(dirO);
                    yT = sin(dirO);
                    C(v,4) = atan(yT,xT);
                    dirExp(v) = dirO;
                else
                    if iPreys(v,3) < collisions then
                        ud(v,:) = zeros(1,2)+uDO;
                        xT = cos(dirO);
                        yT = sin(dirO);
                        C(v,4) = atan(yT,xT);
                    dirExp(v) = dirO;
                    end
                end
            end

            //Orientation-Atraction rules
            if (DetO>0 & DetA>0 & DetR==0 & DetI==0) then
                if (v > Preys) then
                    ud(v,:) = zeros(1,2)+uO;
                    xT = 0.5*cos(dirO)+0.5*cos(dirA);
                    yT = 0.5*sin(dirO)+0.5*sin(dirA);
                    C(v,4) = atan(yT,xT);
                    dirExp(v) = dirO;
                else
                    if iPreys(v,3) < collisions then
                        ud(v,:) = zeros(1,2)+uDO;
                        xT = 0.5*cos(dirO)+0.5*cos(dirA);
                        yT = 0.5*sin(dirO)+0.5*sin(dirA);
                        C(v,4) = atan(yT,xT);
                    dirExp(v) = dirO;
                    end
                end
            end

            //Out of range
            if (DetO==0 & DetR==0 & DetA==0 & DetI==0) then
                if (v > Preys) then
                    ud(v,:) = zeros(1,2)+uO;
                    xT = cos(dirExp(v));
                    yT = sin(dirExp(v));
                    C(v,4) = atan(yT,xT);
                else
                    if iPreys(v,3) < collisions then
                        ud(v,:) = zeros(1,2)+uDO;
                        xT = cos(dirExp(v));
                        yT = sin(dirExp(v));
                        C(v,4) = atan(yT,xT);
                    end
                end
            end

            //ReportsPredators
            if v > Preys then
                ReportPredators(m,v-Preys,1) = C(v,1);
                ReportPredators(m,v-Preys,2) = C(v,2);
                ReportPredators(m,v-Preys,3) = C(v,4);
                ReportPredators(m,v-Preys,4) = C(v,5);
            else
                ReportPreys(m,v,1) = C(v,1);
                ReportPreys(m,v,2) = C(v,2);
                ReportPreys(m,v,3) = C(v,4);
                ReportPreys(m,v,4) = C(v,5);
            end

            Dm(m,v) = C(v,3);
            Dr(v) = max(Dm(:,v))

            Cpast = C(v,:);
            [cs,t] = mov(ud(v,:)',C(v,:)'); //Vehicle movement
            [r, c] = size(cs);
            C(v,:) = cs(:,c)';              //Next initial condition, actual condition

            if (v <= Preys) then
                if iPreys(v,3) >= collisions then
                    C(v,:) = Cpast;
                end
            end

            for i=1:Preys
                iPreys(i,1) = C(i,1)
                iPreys(i,2) = C(i,2)
            end

            if C(v,1) < 0 | C(v,1) > lims | C(v,2) < 0 | C(v,2) > lims then
                C(v,:) = Cpast;
            end

            //This avoids an infinite increment of radians
            if C(v,4) > (2*%pi) then
                C(v,4) = C(v,4)-(2*%pi);
            elseif C(v,4) < 0 then
                C(v,4) = C(v,4)+(2*%pi);
            end
    end

    drawlater();
    clf();
    scatter([0,0],[0,0]);
    d = 0.18;        //Size of arrow
    RadRobot = 0.06; //Radius of robot

    for v=1:(Tv)
        if v <= Preys then
            t1 = C(v,4)-acos(RadRobot/d);
            t2 = C(v,4)+acos(RadRobot/d);
            xf = [C(v,1),C(v,1)+RadRobot*cos(t1),C(v,1)+d*cos(C(v,4)),C(v,1)+RadRobot*cos(t2)];
            yf = [C(v,2),C(v,2)+RadRobot*sin(t1),C(v,2)+d*sin(C(v,4)),C(v,2)+RadRobot*sin(t2)];
            xfpoly(xf,yf,[-1]);
            if iPreys(v,3) < collisions then
                gce().background = 7;
            else
                gce().background = 3;
            end
            xfarc(C(v,1)-RadRobot,C(v,2)+RadRobot,RadRobot*2,RadRobot*2,0,360*64);
            if iPreys(v,3) < collisions then
                gce().background = 7;
            else
                gce().background = 3;
            end
        else
            t1 = C(v,4)-acos(RadRobot/d);
            t2 = C(v,4)+acos(RadRobot/d);
            xf = [C(v,1),C(v,1)+RadRobot*cos(t1),C(v,1)+d*cos(C(v,4)),C(v,1)+RadRobot*cos(t2)];
            yf = [C(v,2),C(v,2)+RadRobot*sin(t1),C(v,2)+d*sin(C(v,4)),C(v,2)+RadRobot*sin(t2)];
            xfpoly(xf,yf,[-1]);
            xfarc(C(v,1)-RadRobot,C(v,2)+RadRobot,RadRobot*2,RadRobot*2,0,360*64);
        end
    end

    a = get("current_axes");
    a.data_bounds = [0,0; lims,lims];
    a.tight_limits = ["on","on"];
    drawnow();
end //of m

cx = sum(ReportPredators(m,:,1))/Predators;
cy = sum(ReportPredators(m,:,2))/Predators;
dy=sqrt(sum((ReportPredators(m,:,2)-cy).^2)/Predators);
dx=sqrt(sum((ReportPredators(m,:,1)-cx).^2)/Predators);
area = (%pi*dy*dx)
set(gca(),"auto_clear","off");
plot(cx,cy,'rx');
xset("color",1);
xarc(cx-(dx*1.3),cy+(dy*1.3),(dx*1.3)*2,(dy*1.3)*2,0,360*64);
//title('k=250, con influencia');
set(gca(),"auto_clear","on");
xlabel('x-axis','Fontsize',3);
ylabel('y-axis','Fontsize',3);
res(1,1)=cx;
res(2,1)=cy;
res(3,1)=dy;
res(4,1)=dx;
res(5,1)=area;
endfunction
