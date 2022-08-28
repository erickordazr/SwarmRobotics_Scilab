function [norm_val]=normalize(val,l_min,l_max)
    norm_val = (val-l_min)/((l_max-l_min)+%eps);
    if norm_val > 1 then
        norm_val = 1;
    elseif norm_val < 0 then
        norm_val = 0;
    end
endfunction

function [val]=denormalize(norm_val,l_min,l_max)
    val = norm_val*(l_max-l_min)+l_min;
    if val > l_max then
        val = l_max;
    elseif val < l_min then
        val = l_min;
    end
endfunction

function [dxdt]=DynamicModel(t,c)
    //c(1,:) x position
    //c(2,:) y position
    //c(3,:) movement
    //c(4,:) orientation
    //c(5,:) speed
    //c(6,:) angular spee

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

    dxdt = [[cos(c(4)) -d*sin(c(4)); sin(c(4)) d*cos(c(4))]*[c(5); c(6)]; 
    c(5); c(6); 
    inv(M + B*inv(Kl)*Ts*inv(A)) * (B*inv(Kl)*Ks*u - (H + F + B*inv(Kl)*inv(A)*[c(5); c(6)]))];
endfunction

function [c,t]= Mov(u,ci)
    ti = 0;              //Initial time
    tf = 1;              //Final time
    tspan = [ti:0.1:tf]; //Time vector
    c = ode(ci,ti,tspan,DynamicModel);
    t = tspan;
endfunction

function [Report,Oend,Oini,st,Tr,Dr]=SimMov2(Ob,Tv,dRU,dOU,dAU)
    st = 0;                         //Iterations
    Report = zeros(st,Tv,4);        //Report
    Tr = zeros(Tv,3);               //Delivery time,Search time,collected objects
    C = zeros(Tv,6);                //vehicle states
    rm = rand(1,1,'normal')*0.01;   //White noise in DC motors (1%)
    uO = 0; uR = 0; uA = 0; uI = 0; //Initial speeds (voltages)
    ud = zeros(Tv,2)+uO;            //Desired speed (voltage)
    lims = 10;                      //Test area limits
    dNest = 4;                      //Maximum distance of influence (nest)4
    dObBox = 3;                   //Maximum distance of influence (Objects box)2.5

    dR = 0.075+dRU;
    dO = 0.075+dOU;
    dA = 0.075+dAU;

    Gripc = zeros(Tv,1); //Open grip
    NestFull = Ob;       //Nest full (End task)

    //Objects box
    centerBox = 0.75;
    limsBox = 0.2;
    ObBox = [centerBox*lims,centerBox*lims];
    ObBoxOn = zeros(Tv,1);

    Objects = zeros(Ob,2);    //Objects location
    if Ob == 0 then
        Object = zeros(1,1)   //Objects vector
        gov = zeros(1,1)      //Objects gripped by vehicle
    else
        Object = zeros(Ob,1); //Objects vector
        gov = zeros(Ob,1)     //Objects gripped by vehicle
    end

    //Random objects position
    for o = 1:Ob
        Obrand1 = denormalize(rand(1,'uniform'),centerBox-(limsBox/2),centerBox+(limsBox/2));
        Obrand2 = denormalize(rand(1,'uniform'),centerBox-(limsBox/2),centerBox+(limsBox/2));
        Objects(o,1) = [lims*Obrand1()];
        Objects(o,2) = [lims*Obrand2()];
        Oini(o,1) = Objects(o,1);       //Initial position of object on X axis
        Oini(o,2) = Objects(o,2);       //Initial position of object on Y axis
        Object(o,1) = 1;
    end

    //Nest
    limsNest = 0.2;
    dot = zeros(1,2)+lims*(limsNest/2); //Nest dot
    Nest = [dot(1,1),dot(1,2)];         //Nest location (Dotted line)
    NestOn = zeros(Tv,1);               //Influence of nest activated for vehicle

    for v = 1:Tv                        //Iinitialize vehicles
        if v == 1 then
            C(v,1) = lims*denormalize(rand(1,'uniform'),0,limsNest); //x position
            C(v,2) = lims*denormalize(rand(1,'uniform'),0,limsNest); //y position
        else
            while %T
                safe = 0;
                C(v,1) = lims*denormalize(rand(1,'uniform'),0,0.2); //x position
                C(v,2) = lims*denormalize(rand(1,'uniform'),0,0.2); //y position
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
        C(v,3) = 0;                                      //movement made
        C(v,4) = denormalize(rand(1,'uniform'),0,2*%pi); //orientation
        C(v,5) = 0;                                      //speed
        C(v,6) = 0;                                      //angular speed
    end

    dirExp = C(:,4);
    explore = zeros(1,Tv);

    /*scf();
    f = get("current_figure");
    f.BackgroundColor = [1 1 1];
    f.figure_position = [10,10];
    f.figure_size=[700,800];*/

    //Finish task when nest is full
    while NestFull <> 0

        st = st +1; //Iteration counter

        for v = 1:Tv

            //Elements detected
            DetR = 0;     //Repulsion
            DetO = 0;     //Orientation
            DetA = 0;     //Attraction
            DetN = 0;     //Nest and influence (Objects)
            DetObBox = 0; //Object Box
            ElemRx = [];   //Elements in repulsion
            ElemOx = [];   //Elements in orientation
            ElemAx = [];   //Elements in attraction
            ElemRy = [];   //Elements in repulsion
            ElemOy = [];   //Elements in orientation
            ElemAy = [];   //Elements in attraction

            //Detection range for repulsion, orientation, attraction, influence, nest and objects
            rangR = 5.75959;
            rangO = 2.44346
            rangA = 0.523599;
            rangI = 0.523599;
            rangN = 2.44346;
            rangOB = 2.44346;

            //Verify each sensor for repulsion of walls
            for i = 1:5
                if i == 1 then
                    dirObs = C(v,4)-3.8397244;
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
                if limitX > lims | limitX < 0 | limitY > lims | limitY < 0 then
                    explore(v) = 0;
                    dirExp(v) = dirObs+(3*%pi/4)+(rand()*%pi/2);
                    if dirExp(v) > (2*%pi) then
                        dirExp(v) = dirExp(v)-(2*%pi);
                    end
                end
            end

            //Resulting direction due object box
            dirOb = atan(ObBox(1,2)-C(v,2),ObBox(1,1)-C(v,1));
            if dirOb < 0 then
                dirOb = dirOb+(2*%pi);
            elseif dirOb > (2*%pi) then
                dirOb = dirOb-(2*%pi);
            end

            //Calculation of influence angles by object box
            BetaOB = dirOb-C(v,4);
            if  BetaOB < 0 then
                BetaOB = BetaOB+(2*%pi);
            end
            GamaOB = C(v,4)-dirOb;
            if  GamaOB < 0 then
                GamaOB = GamaOB+(2*%pi);
            end
            if GamaOB < BetaOB then
                DeltaOB = GamaOB;
            else
                DeltaOB = BetaOB;
            end

            //Calculated distance between the robots and the object zone
            if DeltaOB < rangOB/2 then
                dOB = sqrt(((C(v,1)-ObBox(1,1))^2)+((C(v,2)-ObBox(1,2))^2));
            else
                dOB = %inf;
            end

            [norm_valOB] = normalize(dOB,0,dObBox)  //normalize distance of influence by nest
            [valOB] = denormalize(norm_valOB,1,2.8) //desnormalize distance of influence by nest
            uI = valOB;

            ///////////////////////////////////////////////////////////////////

            for z = 1:Tv
                if v <> z then //it must not be the same

                    //Angle of the individual with respect to other members of the swarm
                    ang = atan(C(z,2)-C(v,2),C(z,1)-C(v,1))
                    if ang < 0 then
                        ang = ang+(2*%pi);
                    elseif ang > (2*%pi) then
                        ang = ang-(2*%pi);
                    end
                    
                    //Calculation of angles of repulsion and attraction with respect to other individuals
                    Beta = ang-C(v,4);
                    if Beta < 0 then
                        Beta = Beta+(2*%pi);
                    end
                    Gama = C(v,4)-ang;
                    if Gama < 0 then
                        Gama = Gama+(2*%pi);
                    end
                    if Gama < Beta then
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
            end //of z
            
            ///////////////////////////////////////////////////////////////////
            
            for o = 1:Ob
                
                //Search object
                if Object(o,1) == 1 & Gripc(v,1) == 0 then
                    
                    //Angle of influence (respect to objects)
                    angI = atan(Objects(o,2)-C(v,2),Objects(o,1)-C(v,1)); 
                    if angI < 0 then
                        angI = angI+(2*%pi);
                    elseif angI > (2*%pi) then
                        angI = angI-(2*%pi);
                    end
                    
                    //Calculation of influence angles (respect to objects)
                    BetaI = angI-C(v,4);
                    if  BetaI < 0 then
                        BetaI = BetaI+(2*%pi);
                    end
                    GamaI = C(v,4)-angI;
                    if  GamaI < 0 then
                        GamaI = GamaI+(2*%pi);
                    end
                    if GamaI < BetaI then
                        DeltaI = GamaI;
                    else
                        DeltaI = BetaI;
                    end
                    
                    //Noise for distance value in influence (standard distribution) of 5%
                    ds1 = 0;
                    ds2 = 0;
                    for i = 1:12
                        ds1 = ds1+rand(1,1,'normal');
                    end
                    for j = 1:6
                        ds2 = ds2+rand(1,1,'normal');
                    end
                    rds= (ds1-ds2)*0.05;

                    //Calculated influence distance (respect to objects)
                    if DeltaI < rangI/2 then
                        di = sqrt((C(v,1)-Objects(o,1))^2+((C(v,2)-Objects(o,2))^2))+rds;
                    else
                        di = %inf;
                    end
                    
                    //Distance between objects and robot
                    Ilim = 0.02;
                    if di <= Ilim then
                        NestOn(v) = 1;
                        gov(o) = v;
                        Tr(v,3) = Tr(v,3)+1;
                        Gripc(v,1) = 1;  //Close grip
                        Object(o,1) = 0; //Object taken by robot
                    end
                    
                end
                
                //Entrega Nest
                if Object(o,1) == 0 & NestOn(v) == 1 then
                    
                    //Angle respect to nest
                    dirNest = atan(Nest(1,2)-C(v,2),Nest(1,1)-C(v,1));
                    if dirNest < 0 then
                        dirNest = dirNest+(2*%pi);
                    elseif dirNest > (2*%pi) then
                        dirNest = dirNest-(2*%pi);
                    end
                    
                    //Calculation of influence angles
                    BetaN = dirNest-C(v,4);
                    if BetaN < 0 then
                        BetaN = BetaN+(2*%pi);
                    end
                    GamaN = C(v,4)-dirNest;
                    if GamaN < 0 then
                        GamaN = GamaN+(2*%pi);
                    end
                    if GamaN < BetaN then
                        DeltaN = GamaN;
                    else
                        DeltaN = BetaN;
                    end
                    
                    //Calculated influence distance (respect to nest)
                    if DeltaN < rangN/2 then
                        dN = sqrt((Nest(1,1)-C(v,1))^2+(Nest(1,2)-C(v,2))^2);
                    else
                        dN = %inf;
                    end
                    
                    //Distance between nest and robot
                    Nlim = 0.05;
                    if dN < Nlim then
                        Objects(o,1) = Nest(1,1);
                        Objects(o,2) = Nest(1,2);
                        dot(1,1) = dot(1,1) + denormalize(rand(1,'uniform'),-0.1,0.1)
                        dot(1,2) = dot(1,2) + denormalize(rand(1,'uniform'),-0.1,0.1)
                        Nest(1,1) = dot(1,1);
                        Nest(1,2) = dot(1,2);
                        NestFull = NestFull - 1;
                        Gripc(v,1) = 0;
                        NestOn(v) = 0;
                    end
                    
                    [norm_val] = normalize(dN,Nlim,dNest) //normalize distance of influence by nest
                    [valN] = denormalize(norm_val,1.2,2.4) //desnormalize distance of influence by nest
                    uI = valN; 
                    
                    if Gripc(gov(o,1),1) == 1 then
                        Objects(o,1) = C(gov(o,1),1);
                        Objects(o,2) = C(gov(o,1),2);
                    end
                end
                
            end //of o
            
            uR = 1.2;
            uO = uI;
            uA = (uI+2.4);
            
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
            
            //Behavior Policies
            ///////////////////////////////////////////////////////////////////
            
            //Repulsion behavior
            if DetR > 0 then
                if NestOn(v) == 0 then
                    if dOB < dObBox then
                        ud(v,:) = zeros(1,2)+uR+rm;
                        xT = 0.8*cos(dirR)+0.2*cos(dirOb);
                        yT = 0.8*sin(dirR)+0.2*sin(dirOb);
                        C(v,4) = atan(yT,xT);
                        //dirExp(v) = dirR;
                        explore(v) = 1;
                    else
                        ud(v,:) = zeros(1,2)+uR+rm;
                        xT = cos(dirR);
                        yT = sin(dirR);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirR;
                    end
                else
                    if dN < dNest then
                        ud(v,:) = zeros(1,2)+uR+rm;
                        xT = 0.8*cos(dirR)+0.2*cos(dirNest);
                        yT = 0.8*sin(dirR)+0.2*sin(dirNest);
                        C(v,4) = atan(yT,xT);
                    else
                        ud(v,:) = zeros(1,2)+uR+rm;
                        xT = cos(dirR);
                        yT = sin(dirR);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirR;
                    end
                end
            end
            
            //Orientation of vehicles
            if DetR == 0 & DetO > 0 & DetA == 0 then
                if NestOn(v) == 0 then
                    if dOB < dObBox then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirOb);
                        yT = sin(dirOb);
                        C(v,4) = atan(yT,xT);
                        explore(v) = 1;
                    else
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirO);
                        yT = sin(dirO);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirO;
                    end
                else
                    if dN < dNest then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirNest);
                        yT = sin(dirNest);
                        C(v,4) = atan(yT,xT);
                    else
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirO);
                        yT = sin(dirO);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirO;
                    end
                end
            end
            
            //Attraction of vehicles
            if DetR == 0 & DetO == 0 & DetA > 0 then
                if NestOn(v) == 0 then
                    if dOB < dObBox then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirOb);
                        yT = sin(dirOb);
                        C(v,4) = atan(yT,xT);
                        explore(v) = 1;
                    else
                        ud(v,:) = zeros(1,2)+uA+rm;
                        xT = cos(dirA);
                        yT = sin(dirA);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirA;
                    end
                else
                    if dN < dNest then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirNest);
                        yT = sin(dirNest);
                        C(v,4) = atan(yT,xT);
                    else
                        ud(v,:) = zeros(1,2)+uA+rm;
                        xT = cos(dirA);
                        yT = sin(dirA);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = dirA;
                    end
                end
            end
            
            //Orientation and Attraction of vehicles
            if DetR == 0 & DetO > 0 & DetA > 0 then
                if NestOn(v) == 0 then
                    if dOB < dObBox then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirOb);
                        yT = sin(dirOb);
                        C(v,4) = atan(yT,xT);
                        explore(v) = 1;
                    else
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = 0.5*cos(dirO)+0.5*cos(dirA);
                        yT = 0.5*cos(dirO)+0.5*sin(dirA);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = atan(sin(dirO)+sin(dirA),cos(dirO)+cos(dirA));
                    end
                else
                    if dN < dNest then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirNest);
                        yT = sin(dirNest);
                        C(v,4) = atan(yT,xT);
                    else
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = 0.5*cos(dirO)+0.5*cos(dirA);
                        yT = 0.5*cos(dirO)+0.5*sin(dirA);
                        C(v,4) = atan(yT,xT);
                        dirExp(v) = atan(sin(dirO)+sin(dirA),cos(dirO)+cos(dirA));
                    end
                end
            end
            
            //Out of range
            if DetR == 0 & DetO == 0 & DetA == 0 then
                if NestOn(v) == 0 then
                    if dOB < dObBox then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirOb);
                        yT = sin(dirOb);
                        C(v,4) = atan(yT,xT);
                        explore(v) = 1;
                    else
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirExp(v));
                        yT = sin(dirExp(v));
                        C(v,4) = atan(yT,xT);
                    end
                else
                    if dN < dNest then
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirNest);
                        yT = sin(dirNest);
                        C(v,4) = atan(yT,xT);
                    else
                        ud(v,:) = zeros(1,2)+uO+rm;
                        xT = cos(dirExp(v));
                        yT = sin(dirExp(v));
                        C(v,4) = atan(yT,xT);
                    end
                end
            end
            
            if explore(v) == 1 & dOB > dObBox & rand(1,'uniform') < 0.1 then
                explore(v) = 0;
                dirExp(v) = dirExp(v)+(3*%pi/4)+(rand()*%pi/2);
                if dirExp(v) > (2*%pi) then
                    dirExp(v) = dirExp(v)-(2*%pi);
                end
            end
            
            ///////////////////////////////////////////////////////////////////
            
            //Reports
            Report(st,v,1) = C(v,1);
            Report(st,v,2) = C(v,2);
            Report(st,v,3) = C(v,3);
            Report(st,v,4) = C(v,4);
            
            Dm(st,v) = C(v,3);
            Dr(v) = max(Dm(:,v))
           
            Cpast = C(v,:);
            [cs,t] = Mov(ud(v,:)',C(v,:)'); //Vehicle movement
            [r, c] = size(cs); 
            C(v,:) = cs(:,c)';              //Next initial condition, actual condition
            
            if C(v,1) < 0 | C(v,1) > lims | C(v,2) < 0 | C(v,2) > lims then
                C(v,:) = Cpast;
            end
            
            //This avoids an infinite increment of radians
            if C(v,4) > (2*%pi) then
                C(v,4) = C(v,4)-(2*%pi);
            elseif C(v,4) < 0 then
                C(v,4) = C(v,4)+(2*%pi);
            end
            
            //Delivery time
            if Gripc(v,1) == 1 then //Close grip
                Tr(v,1) = Tr(v,1) + 1;
            end
           
       end //of v
        
        animation(C,Tv,Objects,Ob,lims,limsNest,centerBox,limsBox,dNest,dObBox);
        
    end //of while
    
    Tr(:,2) = st - Tr(:,1);
    Oend(o,1) = Objects(o,1);
    Oend(o,2) = Objects(o,2);
    
endfunction

function animation(C,Tv,Objects,Ob,lims,limsNest,centerBox,limsBox,dNest,dObBox)
    
    drawlater();
    clf();
    scatter([1,1],[1,1]);
    xrect(0,limsNest*lims,limsNest*lims,limsNest*lims);
    gce().foreground = color("red");
    xarc(lims*(limsNest/2)-dNest,lims*(limsNest/2)+dNest,2*dNest,2*dNest,0,360*64);
    gce().foreground = color("red");
    xrect(lims*(centerBox-(limsBox/2)),lims*(centerBox+(limsBox/2)),lims*limsBox,lims*limsBox);
    gce().foreground = color("blue");
    xarc(lims*centerBox-2.5,lims*centerBox+2.5,2*2.5,2*2.5,0,360*64);
    gce().foreground = color("blue");
    
    d = 0.3;        //Size of arrow
    RadRobot = 0.1; //Radius of robot
    for v = 1:Tv
        t1 = C(v,4)-acos(RadRobot/d);
        t2 = C(v,4)+acos(RadRobot/d);
        xf = [C(v,1),C(v,1)+RadRobot*cos(t1),C(v,1)+d*cos(C(v,4)),C(v,1)+RadRobot*cos(t2)];
        yf = [C(v,2),C(v,2)+RadRobot*sin(t1),C(v,2)+d*sin(C(v,4)),C(v,2)+RadRobot*sin(t2)];
        xfpoly(xf,yf,[-1]);
        xfarc(C(v,1)-RadRobot,C(v,2)+RadRobot,RadRobot*2,RadRobot*2,0,360*64);
    end
    
    RadOb = 0.1; //Radius of object
    for o = 1:Ob
        xfarc(Objects(o,1)-RadOb,Objects(o,2)+RadOb,RadOb*2,RadOb*2,0,360*64);
        gce().background = color("green");
    end
    
    a = get("current_axes");
    a.data_bounds = [0,0; lims,lims];
    a.tight_limits = ["on","on"];
    drawnow();
    
endfunction

function [rst,rDr,rTr1,rTr2,Report]=replicas(Ob,Tv,dRU,dOU,dAU)
    for i=1:3
        [Report,Oend,Oini,st,Tr,Dr]=SimMov2(Ob,Tv,dRU,dOU,dAU)
        Dr = mean(Dr);
        Tr1 = mean(Tr(:,1));
        Tr2 = mean(Tr(:,2)); 
        for j=1:1
            rst(i,j) = st;
            rDr(i,j) = Dr; 
            rTr1(i,j) = Tr1;
            rTr2(i,j) = Tr2;
        end
    end
    prst = mean(rst);
    prDr = mean(Dr); 
    prTr1 =mean(rTr1);
    prTr2 =mean(rTr2);
    Report = [prTr1 prTr2 prst prDr]
endfunction
